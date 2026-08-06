# Chemical Probes Portal -> probedb

Start here. This directory turns one export from
[chemicalprobes.org](https://www.chemicalprobes.org) into a staging directory and
loads it into the database. It is the worked example of §3 of the repository
README, done against a real source rather than the template.

```bash
uv venv --python 3.12 .venv && uv pip install -e . && uv pip install pandas loguru rdkit pytest

python chemicalprobes.org/preprocess.py \
    --json chemicalprobes.org/ChemicalProbesPortal-6_8_2026.json \
    --out  chemicalprobes.org/staging \
    --db   probe.db
```

Reads only the JSON, writes only the staging directory, about four seconds. No
network: the one lookup it needs was done once and is checked in.

## What it changed in your files

**`database/schema.sql`** — 7 tables and 2 views appended at the end, under a
banner comment. **`database/probedb/schema.py`** — their names added to `TABLES`
and `VIEWS`, and one regex anchored (below). Nothing else outside this directory:

```bash
git diff db_prt..chemicalprobes.org -- loader database/probedb/db.py examples tests.py
# prints nothing
```

The tables are **appended, not inserted**: `probedb/schema.py` reads a vocabulary
back out of `schema.sql` with a regex, and a `CHECK` on a column whose name ends
in `type` or `relation` placed above the table it belongs to would be found first
and silently replace that table's vocabulary. The regex is now anchored on a word
boundary, which removes the hazard, but appending stays the habit.

## The export, and what it maps onto

1247 probes, 18 keys each, six of them lists — **34 key paths**. The notebook
asserts that its mapping table and a recursive path walk of the JSON name the
same 34, so a field cannot go missing silently.

| portal | schema | note |
| ------ | ------ | ---- |
| `name`, `InChIkey`, `smiles` | `compound` | 1190 of 1191 SMILES recompute to the stated key |
| `ChEMBL_ID[]` | `chembl` | `\|`-joined in one cell, the loader splits it |
| `canSAR_ID` | `compound_xref` | 1214 values, 1214 distinct |
| `primary_targets[].uniprot_id` | `uniprot` | 656 accessions, all canonical |
| `primary_targets[].name` | `uniprot.hgnc`, reused as `target.name` | **HGNC symbols**, not protein names |
| `primary_targets` | `target.type = 'protein'` | one accession each, so no complex or family |
| `primary_targets[].class` / `subClass` | `target_class` | a family taxonomy, **not** `target.type` |
| `moa` | `bioactivity.moa` | 35 spellings collapse to 30 |
| `inCellValidations[].potency` | `bioactivity_type` | the endpoint label, not a number |
| `*Validations[].potencyValue` | `relation` / `value` / `unit` | free text, parsed |
| `*Validations[].assayDesc` | `assay_description` | |
| the in-vitro / in-cell tier | `assay_type` = `biochemical` / `cell` | the portal's "in vitro" means cell-free |
| `unsuitable`, `pains`, `toxicophore`, ratings, `published_date` | `probe_assessment` | |
| `inVivoValidations[]` | `in_vivo_dose` | a dose is not a potency |
| `PMID[]` | `compound_reference` | named after PubMed, holds DOIs and bare urls |
| `control_compounds[]` | `compound_annotation` | names only, no structures |
| `URL` | `bioactivity_source.xref_id` + `source_xref` | |

**Four names mislead.** `primary_targets[].name` is a gene symbol, not a protein
name. `class` is a taxonomy, not `target.type`. `organism` is the dosed animal,
not `uniprot.species`. `PMID` is a reference list in four formats, three of them
with typo'd hosts.

## What comes out

11 files, every one a table, loaded row for row.

| file | rows | table |
| ---- | ---: | ----- |
| `new_chemicalprobes_compound.tsv` | 1223 | `compound`, and `chembl` from the `\|`-joined ids |
| `new_chemicalprobes_target.tsv` | 644 | `target` |
| `new_chemicalprobes_uniprot.tsv` | 644 | `uniprot`, and `target_uniprot` from `target_key` |
| `new_chemicalprobes_bioactivity.tsv` | 3626 | `bioactivity`, `bioactivity_source`, `bioactivity_group` |
| `new_chemicalprobes_probe_assessment.tsv` | 1223 | `probe_assessment` |
| `new_chemicalprobes_compound_xref.tsv` | 1214 | `compound_xref` |
| `new_chemicalprobes_target_class.tsv` | 1475 | `target_class` |
| `new_chemicalprobes_in_vivo_dose.tsv` | 1450 | `in_vivo_dose` |
| `new_chemicalprobes_compound_reference.tsv` | 1792 | `compound_reference` |
| `new_chemicalprobes_compound_annotation.tsv` | 421 | `compound_annotation` |
| `new_chemicalprobes_rejected_record.tsv` | 47 | `rejected_record` |

**17366 rows**, `PRAGMA foreign_key_check` clean. The first four are the staging
contract and load through `loader/` unchanged; the rest load through
`preprocess.py`, because `loader/` reads the four files of the contract and
nothing else.

## Why each new table exists

| table | why it is not a column somewhere |
| ----- | -------------------------------- |
| `probe_assessment` | one source's verdict per compound. `verdict IN ('recommended','unsuitable')` so a source can say either; without it, 260 of the 261 compounds with no measurement are indistinguishable from a failed download |
| `compound_xref` | the `chembl` treatment generalised: an identifier in another resource, `resource` a CHECK vocabulary so the next one is a value and not a migration |
| `target_class` | a source's taxonomy, keyed on `target_id`. Both levels are multi-valued — 215 of 1475 values would be lost as columns |
| `in_vivo_dose` | a dose has no target and `bioactivity.target_id` is `NOT NULL`; 311 rows could not pick one, and `mg/kg` is not a concentration |
| `compound_reference` | a reading list per compound, not a paper per measurement. `raw` is unrecoverable in 88% of rows |
| `compound_annotation` | the escape hatch, and it should stay small: 421 rows, only what has no column anywhere |
| `rejected_record` | one curator worklist — a number that could not be trusted, or a record that could not become a compound at all |

Views: `probe_flat` (the curator's sheet in one SELECT) and `target_class_flat`.

## Two bugs in `probedb`/`loader` — found, not fixed

Both are one-line changes, both yours to take or leave. Two tests pin the current
behaviour, each with a comment saying to delete it once fixed.

**`database/probedb/db.py:297`** — `unit=unit,` is missing the `or None` every
neighbouring column has, so a blank unit stores as `''` while `loader/load.py:82`
looks for `None`. Every row without a unit is re-inserted on every load:
`3593 → 3742 → 3891`. Fix: `unit=unit or None`.

**`loader/load.py:60,82`** — the duplicate identity omits `assay_description`,
`assay_type` and `cell_line`, so two different assays that agree on the number
collapse: 33 rows, only 6 of them genuine repeats. Discarded: a WT and a C1156Y
mutant at 0.6 nM, an ITC and a TR-FRET number at 30 nM, two NMT orthologues at
2 nM. Until it is fixed, `populate()` inserts the collapsed rows itself, so all
3626 reach the database.

## Reading the numbers

`bioactivity.value` is the number **as the source reported it**, nothing
converted. Three things follow:

- **1856 rows have no `bioactivity_type`** — the portal's in-vitro tier carries no
  endpoint. `explore_db.ipynb` groups by it and pandas drops NaN keys, so those
  rows vanish from that summary; pass `dropna=False`.
- **A range is written as its low end with `relation = '>='`**, the range itself
  in `assay_description`. `>=` was already in `ck_relation`.
- **`bioactivity.source`** names the paper (`DOI:…`, `PMID:…`) when a probe cites
  exactly one, the probe otherwise — the convention `staging/_template` uses.
  1984 of 3626 rows name a paper.

## Still open

- **D8** — 103 ranges and 284 ± errors survive only as text in
  `assay_description`. `value_high` and `value_error` would give them columns.
- **D14** — 97 probes repeat a validation verbatim across their paralogues (148
  copies, 165 rows). That is one experiment against a group, which is what
  `target.type = 'family'` is for; loading them as separate proteins is faithful
  to the export and misleading about independence.
- **23 rows in `rejected_record`** need a curator: a factor of `10^n`, a
  reciprocal rate constant, two numbers behind a `/`.

## Tests

```bash
python -m pytest chemicalprobes.org -q            # 377 tests, ~3 s
python -m pytest chemicalprobes.org -q -m "not slow"
```

Written before the implementation existed, against the docstrings in
`preprocess.py`, and red with `NotImplementedError` until it did.
`test_values.py` is the one worth reading: 36 regression cases, each with the
wrong answer an earlier parser gave.

## Where to read more

| | |
| --- | --- |
| `preprocess.ipynb` | the analysis: every field mapped, checked against a path walk of the JSON |
| `REVIEW_FINDINGS.md` | what five adversarial passes found, and what is still open |
| `handoff/` | a frozen bundle of the **previous** 15-table schema, already distributed. Self-contained; do not expect it to match this branch's `probedb` |
| `preprocess.py` | the pipeline. Read `preprocess()` first, it is 20 lines |
