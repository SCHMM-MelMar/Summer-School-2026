# Chemical Probes Portal -> probedb

What `ChemicalProbesPortal-6_8_2026.json` contains, which of it `database/schema.sql`
can hold today, and what has to be decided before anything is written to
`staging/`.

```
ChemicalProbesPortal-6_8_2026.json   the export, 1.33 MB, 1247 probes
preprocess.ipynb                     the analysis, every count below comes from it
preprocess.py                        end-to-end writer, skeleton only so far
unsuitable.tsv                       the 260 unsuitable probes, keyed on InChIKey
PREPROCESSING.md                     this file
REVIEW_FINDINGS.md                   what three adversarial review passes found
```

Nothing has been written to `staging/`, and **no file in `database/` or `loader/`
has been changed** — that is a decision, not an omission, see below.

## Decisions taken

| | decision |
| --- | -------- |
| **schema** | `database/schema.sql` is **not changed**. D2, D3, D5 and D11 are not taken, so the 11 fields with no column, the 1265 in vivo records, and the range/error numbers do not reach the database. The summary of exactly what that omits is below |
| **loader** | `database/` and `loader/` are **not touched**. The two defects stay documented in `REVIEW_FINDINGS.md`; the 156 rows re-inserted per reload and the 35 dropped rows are a property of loading this source as things stand |
| **the keyless 34** | resolved externally. 10 of 34 resolved through PubChem, all 10 round-tripping through RDKit and none colliding with an existing key. Written with the lookup recorded, so a resolved structure is never mistakable for a portal one |
| **the unsuitables** | their own entry, `unsuitable.tsv`, keyed on `inchikey` so it references `compound` the way `chembl` does. 260 rows, 14 columns |

### unsuitable.tsv

```
inchikey  name  smiles  chembl_id  cansar_id  portal_path  published_date
pains  toxicophore  rating_in_cell  rating_in_organism  rating_count
reference  source_db
```

The unsuitables are the one group where a flat entry loses nothing. All 260 carry
a structure, a canSAR id, a portal path, a date and both structural alerts; none
carries a target, a validation or an in vivo record; 99 carry a ChEMBL id and 4
carry references. Their three rating columns are `0` on every row — the portal
does not rate what it has ruled out — so the entry is complete rather than a
subset. Every `inchikey` in it is also in `compound.tsv`, and none repeats.

Because the schema is unchanged there is no table for it to load into yet; it is
the exact contents of one if that changes.

## Files added

| file | what changed |
| ---- | ------------ |
| `preprocess.ipynb` | was empty. The analysis: flatten, field inventory, a mapping table that checks itself against a path walk of the JSON, per-table checks, the `potencyValue` parser with its coverage and residuals, the structure check against RDKit, and 17 decisions |
| `preprocess.py` | was empty. The skeleton of the writer: every function it needs, its contract, and which decision blocks it. No bodies |
| `PREPROCESSING.md` | this file |
| `REVIEW_FINDINGS.md` | the review record: what was wrong in the first pass, what was fixed, what is still open |

Environment, with `uv`:

```bash
uv venv --python 3.12 .venv
uv pip install pandas loguru nbconvert jupytext ipykernel rdkit
uv pip install -e .
```

RDKit is only needed for the structure check in section 4; the notebook skips
that cell if it is missing.

## What is in the export

One object per probe, 18 keys of which six are lists — **34 key paths in all**.
Flattened:

| frame | rows | from |
| ----- | ---- | ---- |
| probe | 1247 | the probe record |
| target | 1372 | `primary_targets[]` |
| validation | 3551 | `inVitroValidations[]` (1652) and `inCellValidations[]` (1899) under a target |
| invivo | 1265 | `inVivoValidations[]`, probe level |
| chembl | 678 | `ChEMBL_ID[]` |
| reference | 1816 | `PMID[]` |
| control | 418 | `control_compounds[]` |

The notebook asserts that the mapping table and a recursive path walk of the JSON
name the same 34 paths, so a field cannot go missing silently. Of the 34:
**16 map to a column today, 4 need a decision first, 11 have no column at all,
3 are not data.**

Three structural facts worth knowing before reading the mapping:

- `inVitroValidations` **on the probe** is empty in all 1247 records. The
  in-vitro numbers live under the primary target, under the same key name.
- 260 probes list no primary target at all, and they are *exactly* the 260 marked
  `unsuitable = "Yes"`. The portal strips target and validation data from
  compounds it has judged unsuitable as probes. Those two facts are one fact,
  which is why dropping the `unsuitable` flag would leave 260 compounds in the
  database that look like a failed download rather than a verdict. They still
  carry 99 ChEMBL ids and 8 references.
- the nesting levels are not uniform. An in-vitro validation has **no `potency`
  key**, and 7 of the 1265 in vivo records have **no `dose` key**, so
  `record["dose"]` raises on them.

## The mapping

Section 3 of the notebook is the full table with the reasoning per path. The
reasoning matters most where the portal's word and the schema's word are the same
word for different things.

### Maps to a column today

| portal | schema |
| ------ | ------ |
| `name` | `compound.name` |
| `InChIkey` | `compound.inchikey` |
| `smiles` | `compound.smiles` |
| `ChEMBL_ID[]` | `chembl.chembl_id` |
| `primary_targets[].uniprot_id` | `uniprot.uniprot_id` |
| `primary_targets[].name` | `uniprot.hgnc`, and reused as `target.name` |
| `primary_targets` | `target.type = 'protein'` |
| `primary_targets[].moa` | `bioactivity.moa` |
| `inCellValidations[].potency` | `bioactivity.bioactivity_type` |
| `*Validations[].potencyValue` | `bioactivity.relation` / `value` / `unit` |
| `*Validations[].assayDesc` | `bioactivity.assay_description` |
| the in-vitro / in-cell tier | `bioactivity.assay_type` = `biochemical` / `cell` |
| `URL` | `bioactivity_source.xref_id` + `source_xref` |

### Where the names mislead

**`primary_targets[].name` is not `target.name`.** The values are HGNC gene
symbols (`BRD4`, `MDM2`, `GSK3B`), so the symbol belongs in `uniprot.hgnc`.
`target.name` in the template holds a protein name (`Bromodomain-containing
protein 4`), which the export does not carry, so the symbol is reused as the
label. Symbol and accession are in strict 1:1 agreement across all 1372 entries.
Note that `add_target` matches on the accession set and never updates the name,
so **whichever source loads first owns the label**: load the portal before a
source with real protein names and `O60885` reads `BRD4` for good.

**`primary_targets[].class` is not `target.type`.** `type` is a `CHECK`ed
vocabulary about composition: `protein`, `complex`, `ppi`, `family`. `class` is a
protein family taxonomy: `Kinase`, `Epigenetic`, `GPCR`, 56 values with 329
subclasses under them, free text and partly typos (`Programmed Cell Deatch`,
`Transcriptional Factor` next to `Transcription factor`). A different axis. It is
also **not a property of the accession**: 69 of 656 accessions are filed under
more than one class and 126 under more than one subclass, because the value
belongs to the portal's (probe, target) curation row. P51812 is both `Kinase` and
`Protein kinase`; Q13822 carries four subclasses.

**`inVivoValidations[].organism` is not `uniprot.species`.** `species` is the
species of the target protein; `organism` is the animal the probe was dosed in
(Mouse 498, Rat 363, Dog 162). And `dose` is not `concentration` either: `mg/kg`
is a dose per body weight, `mg` and `ug` are amounts.

**`PMID` is not a PubMed ID.** A reference list in four formats: 1098 entries
yield a DOI, 580 a PubMed ID, 138 are bare publisher URLs with no identifier to
lift out, and three have a typo'd host (`www.doi/org/`, `www.di.org`). It hangs
off the probe, not off a measurement. 254 probes list the same paper twice, once
as a DOI and once as a PubMed URL.

**The portal's "in vitro" is the schema's "biochemical".** The portal splits
evidence into in-vitro, in-cell and in-organism, where in-vitro means cell-free.
"In vitro" in the wider literature includes cell assays, so the tier name cannot
be copied across as a value. 53 of the 1652 in-vitro descriptions mention a cell
anyway, so the portal's own tiering is not perfectly clean.

**`potency` is a label, `potencyValue` is the number.** `potency` holds `IC50`,
`DC50`, `Dmax` — that is `bioactivity_type`.

## Does the SMILES agree with the InChIKey?

Both go into `compound` and the InChIKey is the primary key, so a disagreement
would make the row self-contradictory. Recomputed with RDKit for all 1191 probes
carrying both:

- every SMILES parses
- **1190 of 1191 recompute to exactly the given InChIKey**, so the pairing can be
  trusted and there is no reason to recompute keys on load
- one mismatch: `NSC117907` is given `JPOAXWSMFOLMQH-UWWJMHSNSA-N` but its SMILES
  hashes to `JPOAXWSMFOLMQH-UHFFFAOYSA-N`. Same skeleton, and `UHFFFAOYSA` is the
  hash of an empty stereo layer, so the SMILES carries no configuration while the
  key claims one — the geometry is on the quinoid `C=N`, written without it. The
  key is the more specific of the two and it is the primary key, so the SMILES is
  the lossy column here.
- **29 of 1191 SMILES are not in RDKit canonical form.** All 29 are the same
  molecule after re-canonicalising, so this is only how the string was written.
  It matters for one thing: comparing `compound.smiles` between sources as
  strings (D16).

No InChIKey repeats, but four pairs share a 14-character skeleton. Three are
genuine stereoisomer pairs that must stay separate (`Crizotinib` /
`S-crizotinib` and one CHEMBL-named pair are the enantiomer hashes
`GFCCVEGC`/`LBPRGKRZ`; another pair is E/Z). The fourth is not: `Intedanib` and
`Ninetedanib` are two misspellings of nintedanib entered twice, with different
ChEMBL and canSAR ids and identical structures once stereochemistry is ignored
(D17).

## Parsing potencyValue

Of the 3551 records, 3316 hold exactly one number, 110 hold several (up to 8),
118 are blank and 7 are text with no number (`submicromolar`, `low uM`,
`Helicase activity`). **3611 rows come out.**

| | rows |
| --- | --- |
| total | 3611 |
| with an operator | 112 |
| with a range | 103 |
| with a ± error | 284 |
| with an assay concentration lifted out of the text | 16 |
| unit inherited from a sibling fragment | 26 |
| no unit at all | 34 |
| quarantined for hand curation | 23 |

Units after normalising: nM 3242, uM 137, % 133, degC 37, none 34, M-1s-1 8,
fold 6, pM 5, h 4, -log(M) 3, min-1 1, M 1.

This is the **second** version of the parser. An adversarial pass over all 3551
records broke the first one in eight ways, four of them producing numbers wrong
by a factor of 10⁶–10⁹ and seven turning a censored value into an exact potency.
The full before/after table is in the notebook and in `REVIEW_FINDINGS.md`. The
worst was `'14 WT, 2.2 T790M, 1.5 L858R, 0.13 L858R/T790M nM'`, where the `M` of
`T790M` matched before the real `nM` and all four values came out in molar.

What cannot be read reliably is now **quarantined rather than guessed**: a factor
of `10^n` the parser refuses to apply, a reciprocal rate constant, two numbers
behind a `/`, a zero molar potency, a range that reads high to low. Those rows
keep their number and their raw text and carry a reason, so they can be curated
by hand. 23 rows.

## Two defects in the existing loader

Both only show up on data like this, where most rows have no unit and no
operator, and both were found by building the staging files and loading them for
real.

| where | what | cost on this data | fix |
| ----- | ---- | ----------------- | --- |
| `database/probedb/db.py:297` | `unit=unit,` is missing the `or None` that every neighbouring column has, so a blank unit is stored as `''` while `loader/load.py` looks for `None` | **156 rows re-inserted on every reload** — the 124 value-less rows plus the 32 that have a value but no unit | `unit=unit or None` |
| `loader/load.py:60,82` | the duplicate identity is `(inchikey, target_id, source_id, bioactivity_type, relation, value, unit)` — `assay_description`, `assay_type` and `cell_line` are not in it, so an ITC and a TR-FRET number that agree, or a Dmax in two cell lines, collapse | **35 rows dropped, only 6 of them genuine repeats** | add `assay_description` and `cell_line` to the identity |

The template has no unitless rows and no colliding assays, which is why neither
has surfaced before.

## The same measurement under several targets

97 probes repeat a validation record verbatim under 2 or more of their primary
targets: 106 distinct records, **148 redundant copies, 165 bioactivity rows**.
MZ1 is the clearest case — `EC50 = 50 nM`, `Inhibit proliferation of MV4;11
cells` appears under BRD3 and again under BRD2, the two differing only by a
trailing tab. An antiproliferation assay measures the compound, not one
bromodomain.

They land on different `target_id`s, so nothing in the schema or the loader can
catch them, and the database looks like it holds independent per-paralogue data.
88 probes list 3 or more paralogues (BRD2/3/4/BRDT, JAK1/2/3, FGFR1-4,
PARP1/PARP2/TNKS/TNKS2), which is what `target.type = 'family'` exists for (D14).

## What would be written

| file | rows |
| ---- | ---: |
| `compound.tsv` | 1223 — 1213 from the portal plus the 10 resolved externally |
| `target.tsv` | 644 |
| `uniprot.tsv` | 644 |
| `bioactivity.tsv` | 3649 — 23 of them quarantined and held back, the rest loadable |
| `unsuitable.tsv` | 260 |

Of those, **35 are silently dropped by the loader as it stands**, only 6 of them
genuine repeats — D13 was declined, so that stands.

## What we are omitting

With the schema unchanged, **16637 records do not reach the database**:

| what | records | why |
| ---- | ------: | --- |
| in vivo records | 1265 | organism and dose, 1215 with a dose. `bioactivity.target_id` is NOT NULL and a dose is not a concentration (D5) |
| references | 1816 | 1098 DOIs, 580 PubMed ids, 138 raw urls. A bioactivity row takes one `source_id` and the export does not say which paper backs which number (D4) |
| `class` values | 1372 | 56 distinct. Not `target.type`, and not a property of the accession (D3) |
| `subClass` values | 1321 | 329 distinct (D3) |
| `control_compounds` | 418 | names only, no structures, and `compound` is keyed on InChIKey (D2) |
| ratings | 3741 | 3 fields per probe. All 0 for the 260 unsuitable, so they carry information only for the other 987 (D2) |
| `pains` / `toxicophore` | 2494 | 91 and 247 `Yes`. Kept for the 260 unsuitable in `unsuitable.tsv` (D2) |
| `published_date`, `canSAR_ID`, `URL` | 3741 | kept for the 260 unsuitable in `unsuitable.tsv` (D2) |
| ± errors | 284 | no numeric column on `bioactivity` (D8) |
| range highs | 103 | no numeric column on `bioactivity` (D8) |
| probes still without a key | 24 | with 31 target entries and 87 measurements under them (D1) |
| rows the loader drops | 35 | only 6 are genuine repeats; D13 declined (D13) |
| rows held for curation | 23 | written to `quarantine.tsv`, not to `bioactivity.tsv` (D15) |

Read the other way: what the four staging files plus `unsuitable.tsv` carry is
every compound, every target, every measurement the portal reports, and the full
portal record for the 260 compounds it has ruled out. What they do not carry is
the portal's *opinion* of the other 987 — its scores, its alerts, its dates and
its reading list — the animal data, and the protein family taxonomy.

Two things about writing the files, both confirmed by writing them and reading
them back through `loader.validate`:

- the free text contains tabs, newlines and quotes (5 tabs, 271 line breaks and
  10 quotes survive stripping in `assayDesc`), so the writer must be
  `csv.writer(delimiter='\t')` and not a `'\t'.join`, which corrupts 130+ rows
- a pandas `NaN` is truthy, so `value or ""` writes the literal string `'nan'`

`validate` will report one line per row without an operator and one per row
without a unit — about 3500 and 156 here. None is a hard error: the portal writes
an operator on 112 of 3611 numbers, and the 124 value-less rows have no unit
either.

## Decisions for review

Nothing is written until these are settled. Full versions with counts are section
14 of the notebook.

| | issue | proposal |
| --- | ----- | -------- |
| **D1** | 34 probes have no InChIKey and no SMILES, so the key that is the PK and every FK cannot be computed | hold them out of this release and list them by name. Costs 34 compounds, 41 target entries, 119 validations, 58 in vivo records, 35 references, 1 ChEMBL id, 1 canSAR id |
| **D2** | 9 probe-level fields plus two lists have no column: 3 ratings, `unsuitable`, `pains`, `toxicophore`, `published_date`, `canSAR_ID`, `URL`, control names, references | one additive table rather than nine columns: `compound_annotation(inchikey, source_db, property, ordinal, value)`, primary key on all four key columns so a reload is idempotent and a list can hold several values. `property`, not `key`, which is reserved in some dialects. ~13100 rows |
| **D3** | `class` / `subClass` have no column, are not `target.type`, and are not a property of the accession either | a parallel `target_annotation` keyed on **`target_id`**, not `uniprot_id`. Keying on the accession, or adding columns to `uniprot`, would let one probe record's typo win globally — and `uniprot` is what `target_flat` exposes. ~1500 rows |
| **D4** | `PMID` is a probe-level reference list, and a bioactivity row takes exactly one `source_id` | `source_db='Chemical Probes Portal'`, `source`=the probe, `xref_id='https://www.chemicalprobes.org/'`, `source_xref`=**the full path**, not the last segment: 260 probes live under `/unsuitables/`. Keep the references as annotations rather than guessing which paper backs which number |
| **D5** | 1265 in vivo records have no target, `bioactivity.target_id` is NOT NULL, and 263 sit under multi-target probes | a small `in_vivo(id, inchikey, organism, dose_value, dose_unit, route, dose_raw, source_id)` table, one row per dose, with a CHECK on `route`. Reject `dose → concentration`: `mg/kg` is not a concentration |
| **D6** | 1785 records have no endpoint — the in-vitro tier has no `potency` key at all, and 133 in-cell records have it null or `Not done` | `bioactivity_type = NULL` and keep the description; optionally backfill the 65 in-vitro records that name it inside `assayDesc`, flagged as derived |
| **D7** | `assay_type` for the 1652 in-vitro rows | `biochemical`, refined to `binding` where `assayDesc` names SPR, ITC, BROMOscan, DSF, MST or a radioligand assay — 379 rows |
| **D8** | 103 ranges and 284 ± errors have no column. Writing the low end with `relation='~'` would claim "approximately 22 nM" about a 22–166 nM range, a stronger claim than the data makes | `>=` is **already** in `ck_relation`, so censoring the low end needs no schema change. Otherwise add `value_high` and `value_error` — two nullable columns that give both a home. Never a computed midpoint |
| **D9** | 110 records hold 2–8 measurements each | one row per number, the fragment text appended to `assay_description` so the qualifying label travels with the number, and the endpoint taken **per fragment** where the record names several. Without that, all 8 rows of one record carry `'DC50, Dmax, IC50'` |
| **D10** | `cell_line` | leave empty. Line names sit inside free text inconsistently and a wrong line is worse than none; Cellosaurus matching is a later enrichment step |
| **D11** | 84 `assayDesc` values exceed the declared `VARCHAR(255)`, the longest 1276 (appending the fragment label per D9 pushes none of the rest over) | widen to `TEXT`. SQLite does not enforce it, but the DDL uses `SERIAL`, so it is written for a database that would |
| **D12** | `uniprot.species` and `entrez_gene` | leave empty. Neither is in the export and the portal is not exclusively human, so `Homo sapiens` would be a guess |
| **D13** | the two loader defects above | fix both before loading anything |
| **D14** | the same measurement under several targets | load as written, one `protein` target each, and record the replication so it is visible. Promoting the replicated records to a `family` target is an inference the export does not make |
| **D15** | quarantined values | write them to their own file with the reason and curate by hand. Add two load-time assertions that would have caught the whole class unaided: reject `unit='M'` with `value > 1`, and reject `value = 0` |
| **D16** | 29 SMILES are not canonical | canonicalise on write with RDKit, or store as given. The InChIKey is the join key either way, so this only buys string comparability |
| **D17** | `Intedanib` and `Ninetedanib` are one compound under two misspellings | load both, since they are two portal records with two InChIKeys, and record the collision. Merging is a curation decision about the source. The other three skeleton pairs are genuine stereoisomers and must stay apart |

**Blocking the writer:** D13 (fix the loader first), then D2, D3, D5 and D11,
which need `database/schema.sql` extended. D1, D8, D9 and D15 change the row
counts above. D6, D7, D10, D12, D16 and D17 only change what sits in a column.

If D2, D3 and D5 are adopted, three things follow: append the tables at the
**end** of `schema.sql`, add them to `TABLES` in `probedb/schema.py` or
`db.counts()` will not see them, and note that `vocabulary()` searches
`schema.sql` with an unanchored regex — a new `CHECK` on a column whose name ends
in `type` or `relation` will hijack an existing vocabulary. Verified DDL for all
three is in `REVIEW_FINDINGS.md`.
