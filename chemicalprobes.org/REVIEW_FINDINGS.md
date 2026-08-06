# Review record

Five independent adversarial passes, each required to reproduce every claim by
running code. Passes 1-3 attacked the analysis and the pipeline; passes 4-5
attacked the schema this branch adds.

## 5. The schema, after it was first written

Seven tables were added quickly and two of them did not survive review.

| was | is | why |
| --- | -- | --- |
| `unsuitable` (260) | `probe_assessment` (1223) | **100% redundant** — every fact was already in `compound_annotation`, and the name was a value masquerading as a table, so a source could never say 'recommended'. Its `source_id` was a fiction: 259 `bioactivity_source` rows with zero measurements |
| `compound_annotation` (11428) | `compound_annotation` (421) + typed columns | 96% of it was a dense 1223x9 rectangle in key/value form. `MAX(value)` on `canSAR_ID` returned **999076** instead of 7447599, because the values are text |
| `target_annotation` | `target_class` | same rows; the README's stated rationale was **false** (`target_id` and `uniprot_id` are a perfect bijection here). The real defence is `ordinal`: 215 of 1475 values are multi-valued |
| `in_vivo` | `in_vivo_dose` | rename only; the shape was right |
| `quarantine` (23) + `skipped_compound` (24) | `rejected_record` (47) | `quarantine` duplicated 9 of `bioactivity`'s 12 columns and had already drifted from it: `TEXT` vs `VARCHAR(255)`, nullable vs `NOT NULL`, and **no `ck_relation` at all**. `skipped_compound` described nothing the database contains |
| — | `compound_xref` | canSAR ids behave exactly like ChEMBL ids: 1214 values, 1214 distinct |

Also fixed: provenance was three conventions in one commit (`source_id`,
`source_db` as text, and nothing at all) and is now `source_id` everywhere;
`bioactivity_source` rows contributing nothing went from **589 to 1**; and
`vocabulary()`'s regex is anchored, so a `CHECK` can no longer shadow another
table's vocabulary from anywhere in the file.

Net: 26784 rows to **17366**, with every count preserved.

---

Passes 1-3, over the first version of this analysis,
each required to reproduce every claim by running code against the 3551 records
rather than reading the notebook. What they found, what was fixed, and what is
still open.

- **pass 1** attacked the `potencyValue` parser for silently wrong numbers
- **pass 2** checked that every key of the export is accounted for, and
  recomputed every quantitative claim
- **pass 3** built the staging files and loaded them into a real SQLite database
  through `loader/`, then judged the proposed schema extensions

Counts below were re-verified independently before being accepted.

## 1. Parser defects, all fixed

Twelve places where the first parser produced a wrong number, unit or operator.
Four were wrong by a factor of 10⁶–10⁹ and fifteen turned a censored or exact
operator into NULL or the wrong sign, which is the worst failure mode for a
potency table.

| raw | was | is | records |
| --- | --- | -- | ------: |
| `'14 WT, 2.2 T790M, 1.5 L858R, 0.13 L858R/T790M nM'` | four values in **M** — the `M` of `T790M` matched before the real `nM`, then propagated to the unitless siblings | four values in nM | 1 (4 rows) |
| `'85n M'` | 85 **M** | 85 nM | 1 |
| `'0. 69 nM'` | **0.0** nM | 0.69 nM | 1 |
| `'6,5 nM'` | **two** rows, 6 and 5 | one row, 6.5 | 1 |
| `'up to 10 µM'`, `'below 1 µM'`, `'about 50 nM'` | relation NULL, censored read as exact | `<`, `<`, `~` | 7 |
| `'IC50 = 493 ± 101 nM'` | relation NULL — the `=` branch of the regex was unreachable | `=` | 8 |
| `'5-550nM'` | upper bound lost | 5–550 nM | 1 |
| `'5 nM to 13 nM; 260 nM and 855 nM'` | 2 of 4 numbers | all 4 | 1 |
| `'760nM and 1000nM'`, `'8.5 and 4.7 nM'`, `'1 or 41 nM'`, `'7.8 nM (human) 2.4 nM (Mouse)'` | second measurement dropped | both kept | 4 |
| `'~50 at 3 uM'` | 50 **uM** — the unit came from the assay concentration | 50 %, concentration 3 uM | 1 |
| `'9900 ± 1800 M–1 s–1'` (en dash, middot) | unit **M** | `M-1s-1`, quarantined | 5 |
| `'63 nM (DC50); 90.8% (Dmax); 52 nM (EC50)'` | all rows labelled `'DC50, Dmax, EC50'` | the endpoint the fragment names, or NULL — never the joined label | 15 (49 rows, 10 of them named) |

Two latent traps were closed by the same guard: a cell line named `(K)` or a
species label `(H)` in a record with no other unit used to produce a temperature
and an hour.

**Confirmed not broken**, so the guards stay: the label-digit lookbehind never
blocks a real value (`BD1`, `MEC1.5`, `H3.3`, `MOLM-13`, `HCT-116` all handled);
no unicode-whitespace failure (NBSP, thin space); no µ/μ/u confusion across the
142 records carrying a micromolar spelling; the `;`-inside-`MV4;11` guard is
needed exactly once; thousands separators are correct.

**What is quarantined instead of guessed** (23 rows): a factor of `10^n`, a
reciprocal rate constant, two numbers behind a `/`, a zero molar potency, a range
that reads high to low, a p-scale endpoint paired with a molar unit, a duration.
Each carries its reason.

## 2. Wrong numbers in the first write-up, all corrected

Every number a code cell printed reproduced exactly. Every error was in
hand-typed prose that had never been recomputed against the frames — which is why
the mapping and the homeless table are now **derived from one list and checked
against a path walk of the JSON**, and no longer typed twice.

| claim | was | is |
| ----- | --- | -- |
| multi-value records | 107 → 163, then 288 | 110 → 295 rows, up to 8 from one record |
| endpoint recoverable from `assayDesc` | 70 | 65 |
| records with no endpoint | 1652 (in-vitro only) | 1785 — 133 in-cell records have `potency` null or `Not done` |
| binding keyword matches | ~260 | 379 |
| ranges | 95 | 103 |
| records holding exactly one number | 3319 | 3316 |
| accessions with >1 subclass | 130 | 126 |
| `assayDesc` pushed over 255 by appending the fragment label | 80 | 0 |
| units census | nM 3235 | nM 3242, uM 137, % 133, degC 37 |
| probes citing one paper as both DOI and PubMed | 322 | 254 |
| in-vitro descriptions mentioning a cell | 55 | 53 |
| `bioactivity_source` rows | ~950 | 1213, one per loadable probe |
| rows re-inserted per reload | 158 | 149 |
| rows the loader drops | 32, 26 distinct | 35, 29 distinct |
| `db.py` line of the unit bug | 301 | 297 |
| subclasses | 330 | 329 — the blank was being counted as a value |
| `compound_annotation` fields | 8 | 9 |
| unsuitable probes contribute "a compound row and nothing else" | — | also 99 ChEMBL ids and 8 references |
| D1 cost | 34/41/119/58 | also 35 references, 1 ChEMBL id, 1 canSAR id |
| `bioactivity.tsv` size | 3488 | 3616 rows in the file |

Three gaps in coverage, all closed:

- `probes[].inVitroValidations` (the probe-level list, empty in all 1247 records)
  appeared in neither the mapping nor the homeless table
- `PMID`, `dose` and `URL` appeared in **both** with contradictory verdicts
- `± error` and `range high` were listed as if they were keys of the export; they
  are values the parser derives

## 3. Loader and schema, verified by loading

The plan loads with **zero hard validation errors and zero constraint
violations** — `loader.load(strict=True)` on a generated staging directory
returned `compounds: 1213, targets: 640` with `PRAGMA foreign_key_check` empty.
Nothing in `schema.sql` has to change to get the data *in*. Three things have to
change to get it in *correctly*. (The bioactivity count in that run was the review
pass's own writer, not the projection in the notebook; the reconciled figure is
3616 written, 35 dropped.)

### Blockers

**`database/probedb/db.py:297`** — `unit=unit,` is missing the `or None` that
every neighbouring column has, so a blank unit is stored as `''` while
`loader/load.py:82` builds its duplicate key with `row.get("unit") or None`. The
key never matches what comes back from the database, so **every row without a
unit is re-inserted on every reload** — 149 of them here: the 124 value-less rows
plus the 31 that carry a value and no unit. After the one-word fix they store as
NULL and a second load re-inserts nothing.

(The review pass reported this at line 301 and as 158 rows, against a staging
directory it generated itself. Re-checked against the file and against the row
counts in the notebook: line **297**, **149** rows.)

**`loader/load.py:60,82`** — the duplicate identity omits `assay_description`,
`assay_type` and `cell_line`, so **35 rows are dropped of which only 6 are
genuine repeats**, losing 29 distinct measurements. Discarded measurements include an ITC and a TR-FRET number
that agree, a Dmax measured in SU-DHL-1 and again in NCI-H2228, a WT and a C1156Y
mutant at the same potency, and two rows that differ only in tier. Adding
`assay_description` and `cell_line` to the identity brings the drop to 6, which
is the right answer; finer source granularity does not help (per-probe sources
still drop the same rows, because `source` is 1:1 with `inchikey` and `inchikey`
is already in the key).

### Confirmed enforced, so the mapping is safe

`ck_relation` rejects `'<<'`, `'≤'`, `''`; `ck_target_type` rejects `'Kinase'`,
`'PROTEIN'`, `''`, which is the hard proof that `class` can never go in that
column; the composite FK to `bioactivity_group(inchikey, target_id, moa)` is
enforced on all three columns and is byte-exact on `moa`. **SQLite enforces no
VARCHAR length at all** — 1276 characters were stored in a `VARCHAR(255)`, and a
60-character InChIKey was accepted into `VARCHAR(27)`. So D11 is cosmetic here and
real on any other engine.

Also verified: **0 of 1331** (compound, accession) pairs carry more than one
`moa`, so the `bioactivity_group` key never splits spuriously.

### Latent hazards in the extension proposals

- **`vocabulary()` uses an unanchored regex** (`probedb/schema.py:23`) with
  `re.search`, taking the first match in the file. Adding a table with
  `annotation_type IN ('probe','target','assay')` **before** `CREATE TABLE
  target` made `vocabulary('type')` return `['assay','probe','target']`, after
  which `validate()` rejected `type='protein'` on a perfectly good directory.
  Append new tables at the end, and consider anchoring the pattern.
- **`schema.TABLES` is hardcoded**, so `db.counts()` silently omits any new table.
- `target.name` is load-order dependent: `add_target` matches on the accession
  set and never updates the name.
- 260 probe URLs have two path segments (`/unsuitables/jib-04`), so
  `rsplit('/')[-1]` builds a 404.

### Verified DDL for D2, D3, D5

Loaded the full portal data against this, appended to the end of `schema.sql`;
`vocabulary()`, `validate()` and `target_flat` all unchanged.

```sql
-- portal/source metadata that is a fact about a compound but has no column.
-- ordinal lets one property hold a list (control compounds, references).
CREATE TABLE compound_annotation (
    inchikey      VARCHAR(27) NOT NULL,
    source_db     VARCHAR(255) NOT NULL,
    property      VARCHAR(100) NOT NULL,
    ordinal       INTEGER NOT NULL DEFAULT 0,
    value         TEXT,
    CONSTRAINT pk_compound_annotation
        PRIMARY KEY (inchikey, source_db, property, ordinal),
    CONSTRAINT fk_ca_compound
        FOREIGN KEY (inchikey) REFERENCES compound (inchikey)
);

-- the same, for facts about a target that are not its composition.
-- keyed on target_id, not uniprot_id: 69 accessions carry more than one
-- portal class and 130 more than one subClass, so it is not a property
-- of the protein.
CREATE TABLE target_annotation (
    target_id     INTEGER NOT NULL,
    source_db     VARCHAR(255) NOT NULL,
    property      VARCHAR(100) NOT NULL,
    ordinal       INTEGER NOT NULL DEFAULT 0,
    value         TEXT,
    CONSTRAINT pk_target_annotation
        PRIMARY KEY (target_id, source_db, property, ordinal),
    CONSTRAINT fk_ta_target
        FOREIGN KEY (target_id) REFERENCES target (target_id)
);

-- a dose given to an animal is not a potency and has no target, so it
-- cannot be a bioactivity row. one row per dose, like a bioactivity
-- fragment; dose_raw keeps the string the dose was read out of.
CREATE TABLE in_vivo (
    id            SERIAL PRIMARY KEY,
    inchikey      VARCHAR(27) NOT NULL,
    organism      VARCHAR(100),
    dose_value    NUMERIC,
    dose_unit     VARCHAR(50),
    route         VARCHAR(20),
    dose_raw      VARCHAR(255),
    source_id     INTEGER,
    source_xref   VARCHAR(100),
    CONSTRAINT ck_in_vivo_route
        CHECK (route IS NULL OR route IN ('IV', 'PO', 'IP', 'SC', 'topical', 'gavage')),
    CONSTRAINT fk_iv_compound
        FOREIGN KEY (inchikey) REFERENCES compound (inchikey),
    CONSTRAINT fk_iv_source
        FOREIGN KEY (source_id) REFERENCES bioactivity_source (source_id)
);
CREATE INDEX idx_in_vivo_compound ON in_vivo (inchikey);
```

Rows these would take: `compound_annotation` ~13100, `target_annotation` ~1500,
`in_vivo` 1207 records excluding the 34 keyless probes, which become 1450 rows
because 193 dose strings hold more than one dose.

## 4. Not taken

`database/schema.sql`, `database/probedb/` and `loader/` are unchanged by
decision. That means the two blockers above are live properties of loading this
source: a reload re-inserts the 149 rows without a unit, and 35 rows are dropped
of which 29 are distinct measurements. Both are one-line fixes whenever they are
wanted, and the counts here are what they cost until then.

## 5. Still open

- **D8** has no fix that is both lossless and free. 103 ranges and 284 ± errors
  have nowhere to go. Adding `value_high` and `value_error` to `bioactivity` is
  two nullable columns and solves it outright; the alternative loses the spread.
  Note that `>=` and `<=` are already in `ck_relation`, so censoring a range's
  low end needs no schema change at all — only the two extra columns would.

- **the multi-endpoint fix is partial.** Assigning the endpoint per fragment
  stops all 8 rows of a record carrying `'DC50, Dmax, IC50'`, but only 10 of the
  49 rows in those 15 records name their endpoint inside the fragment. The other
  39 come out NULL, which is honest and still not the joined label.
- **D14** is a judgement, not a defect. 148 replicated copies across a probe's
  targets are one experiment each; loading them as written is faithful to the
  export and misleading about independence.
- one record still yields a spurious row: `'15 nM; (32 - 46 nM, 0 - 8 h)'`, twice,
  where a parenthesised time window reads as a measurement. It is caught by the
  quarantine rule on duration units rather than by the splitter.
- `db.compound_key()` resolves a lookup against `chembl.chembl_id` before
  `compound.name`, and 18 probes are *named* by a bare ChEMBL id. A lookup by one
  of those names resolves through the wrong table if the id belongs to another
  compound. Not caused by this import, but it will be reachable once it loads.
- 10 probe names carry unstripped whitespace (`'CHEMBL3092538 '`,
  `' (R)-Zinc-3573'`). The writer must strip `name`, which the notebook's
  flattening does not currently do.
