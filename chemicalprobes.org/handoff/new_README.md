# Handoff: Chemical Probes Portal -> probedb

Everything in this directory is prefixed `new_` and is meant to be handed over as
one bundle. Nothing here has to be run: the data is already in table shape.

```
new_schema.sql        database/schema.sql with 7 tables appended, nothing else changed
new_schema.py         database/probedb/schema.py with their 7 names in TABLES
new_<table>.tsv       the built database, one file per table, 15 files, 26784 rows
new_load.py           builds a database from the above, no dependencies
new_README.md         this file
```

Each `new_<table>.tsv` has **exactly the columns its table declares, in order**,
including the surrogate ids, so it loads with no resolution step and nothing to
unpack:

```bash
python3 new_load.py probe.db      # 26784 rows, foreign keys checked
```

Plain Python 3, nothing installed, about two seconds. Use it rather than
sqlite3's own `.import`: `.import` writes an empty cell as the empty string
rather than NULL, and the first bioactivity row then fails
`CHECK (relation IS NULL OR relation IN (...))`, because `''` is neither. I tried
that first and it loaded 0 of 3626 bioactivity rows while reporting success for
everything else.

That is a different shape from the staging files in `../staging/`. Those are what
a *source* hands over — they name a target by its accession and a source by its
name, and `compound.tsv` alone feeds two tables. These are what the *database*
holds. Both describe the same 26784 rows.

The rest of this file is the branch write-up, unchanged.

---

# Chemical Probes Portal -> probedb

Start here. This directory turns one export from
[chemicalprobes.org](https://www.chemicalprobes.org) into a staging directory
and loads it into the database. It is the worked example of §3 of the repository
README, done against a real source rather than the template.

```bash
uv venv --python 3.12 .venv && uv pip install -e . && uv pip install pandas loguru rdkit pytest

python chemicalprobes.org/preprocess.py \
    --json chemicalprobes.org/ChemicalProbesPortal-6_8_2026.json \
    --out  staging/chemicalprobes.org \
    --db   probe.db
```

That is the whole thing. It reads only the JSON, writes only the staging
directory, and takes about four seconds. No network: the one lookup it needs was
done once and is checked in.

## What it changed in your files

**Two files, 148 lines added, nothing removed or edited.**

| file | change |
| ---- | ------ |
| `database/schema.sql` | 7 tables appended at the end, under a banner comment saying so |
| `database/probedb/schema.py` | their 7 names added to `TABLES`, so `db.counts()` sees them |

Nothing else outside this directory is touched. `loader/`, `probedb/db.py`,
`tests.py`, `examples/` and the existing 8 tables are byte-identical to `db_prt`:

```bash
git diff db_prt..chemicalprobes.org -- loader database/probedb/db.py examples tests.py
# prints nothing
```

The 7 tables are **appended, never inserted**, and the banner in `schema.sql`
says why: `probedb/schema.py:23` finds a vocabulary with an unanchored regex, so
a `CHECK` on a column whose name ends in `type` or `relation`, placed above the
table it belongs to, would be found first and silently replace that table's
vocabulary. Verified after the change: `vocabulary("type")` and
`vocabulary("relation")` still return exactly what they did before.

## What comes out

11 files, every one of them a table, loaded row for row.

| file | rows | table |
| ---- | ---: | ----- |
| `chemicalprobes_compound.tsv` | 1223 | `compound`, and `chembl` from the `\|`-joined ids |
| `chemicalprobes_target.tsv` | 644 | `target` |
| `chemicalprobes_uniprot.tsv` | 644 | `uniprot`, and `target_uniprot` from `target_key` |
| `chemicalprobes_bioactivity.tsv` | 3626 | `bioactivity`, `bioactivity_source`, `bioactivity_group` |
| `chemicalprobes_unsuitable.tsv` | 260 | `unsuitable` |
| `chemicalprobes_compound_annotation.tsv` | 11428 | `compound_annotation` |
| `chemicalprobes_target_annotation.tsv` | 1475 | `target_annotation` |
| `chemicalprobes_in_vivo.tsv` | 1450 | `in_vivo` |
| `chemicalprobes_reference.tsv` | 1792 | `compound_reference` |
| `chemicalprobes_quarantine.tsv` | 23 | `quarantine` |
| `chemicalprobes_skipped_compound.tsv` | 24 | `skipped_compound` |

26784 rows, `PRAGMA foreign_key_check` clean. The first four are the staging
contract and load through `loader/` unchanged; the rest load through
`preprocess.py`, because `loader/` reads the four files of the contract and
nothing else and I did not want to change that for one source.

The generated copies are checked in under `chemicalprobes.org/staging/`, so you
can read the output without running anything.

## Why each new table exists

Each one holds something the export carries and the four contract files cannot.

| table | rows | why it is not a column somewhere |
| ----- | ---: | ------------------------------- |
| `unsuitable` | 260 | the portal publishes compounds it has **ruled out**, with no target and no measurement. Without this they are 260 of the 261 compounds in the database with nothing attached and nothing saying why — indistinguishable from a failed download |
| `compound_annotation` | 11428 | scores, alerts, dates, external ids, control-compound names. One property per row so the next source adds one without a migration |
| `target_annotation` | 1475 | the portal's protein-family taxonomy. Keyed on `target_id`, **not** on an accession: 69 accessions carry more than one class, so a column on `uniprot` would let one record's typo win globally |
| `in_vivo` | 1450 | a dose is not a potency — no endpoint, nothing to compare, and no target, so it cannot be a `bioactivity` row. `mg/kg` is also not a concentration |
| `compound_reference` | 1792 | the portal cites a reading list per compound, not a paper per measurement, and a `bioactivity` row takes exactly one `source_id` |
| `quarantine` | 23 | numbers the parser could read but not trust. Kept **out** of `bioactivity` deliberately, and kept rather than dropped so a curator can see them |
| `skipped_compound` | 24 | records that could not become a compound at all. The only table with no foreign key, because its rows are about absence |

Every one carries its foreign keys: `inchikey` onto `compound`, `target_id` onto
`target`, `source_id` onto `bioactivity_source`. A staging file names its source
with `source_db`/`source` and its target with the accession, exactly the way
`bioactivity.tsv` does, and the loading step resolves both.

## Two bugs in `probedb`/`loader` — found, not fixed

Both were found by loading this data and both are one-line changes. I left them
alone deliberately; they are yours to take or leave. Two tests **pin** the
current behaviour so it cannot change unnoticed, each with a comment saying to
delete it once fixed.

**1. `database/probedb/db.py:297` — a blank unit is stored as `''`, looked up as `NULL`.**

```python
relation=relation or None,
value=None if value in (None, "") else float(value),
unit=unit,                        # every neighbour has `or None`; this one does not
assay_type=assay_type or None,
```

`loader/load.py:82` builds its duplicate key with `row.get("unit") or None`, so
the key never matches what comes back. Every row without a unit is re-inserted on
every load — 149 of them here, growing by 149 each time:

```
load 1: inserted 3593   table 3593
load 2: inserted  149   table 3742
load 3: inserted  149   table 3891
```

The README promises "loading the same directory twice does not duplicate
anything". It holds for the template, which has no unitless rows. **Fix:
`unit=unit or None`.**

**2. `loader/load.py:60,82` — the duplicate identity is too narrow.**

It is `(inchikey, target_id, source_id, bioactivity_type, relation, value, unit)`.
`assay_description`, `assay_type` and `cell_line` are not in it, so two different
assays that agree on the number collapse into one — 33 rows here, of which only 6
are genuine repeats:

| discarded | the two assays |
| --------- | -------------- |
| Q9UM73, 0.6 nM | `0.6 nM WT` vs `0.6 nM C1156Y` — a mutant |
| Q92793, 30 nM | `Isothermal calorimetry` vs `TR-FRET` |
| O60551, 2 nM | `2 nM TbNMT` vs `2 nM LmNMT` — two orthologues |

The earlier row wins, which is an accident of file order, not a judgement. Adding
`assay_description` and `cell_line` to both the `SELECT` and the key tuple brings
the drop to exactly the 6 genuine repeats. Until then `preprocess.py` inserts the
collapsed rows itself after `load()`, so all 3626 reach the database.

## Reading the numbers

`bioactivity.value` is the number **as the source reported it**, nothing
converted, which is the schema's own doctrine. Three things follow that are worth
knowing before you query:

- **1856 rows have no `bioactivity_type`.** The portal's in-vitro tier carries no
  endpoint at all. `explore_db.ipynb` groups by that column and pandas drops NaN
  keys, so those rows vanish from that summary — pass `dropna=False`.
- **A range is written as its low end with `relation = '>='`**, and the range
  survives in `assay_description`. `>=` was already in `ck_relation`; no schema
  change was needed for it.
- **`bioactivity.source`** names the paper (`DOI:…`, `PMID:…`) when a probe cites
  exactly one, and the probe otherwise — the convention `staging/_template`
  already uses for Probes & Drugs. 1984 of 3626 rows name a paper.

## Tests

```bash
python -m pytest chemicalprobes.org -q            # 388 tests, ~3 s
python -m pytest chemicalprobes.org -q -m "not slow"   # skips the ones that read the export
```

Written before the implementation existed, against the docstrings in
`preprocess.py`, and red with `NotImplementedError` until it did.
`test_values.py` is the one worth reading: 36 regression cases, each with the
wrong answer an earlier version of the parser gave.

## Where to read more

| | |
| --- | --- |
| `preprocess.ipynb` | the analysis: every field of the export mapped onto the schema, checked against a path walk of the JSON so a field cannot go missing silently |
| `PREPROCESSING.md` | the decisions, with the counts behind each |
| `REVIEW_FINDINGS.md` | what three adversarial review passes found and what is still open |
| `preprocess.py` | the pipeline. Read `preprocess()` first, it is 20 lines |
