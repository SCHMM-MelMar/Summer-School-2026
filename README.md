# DB

A small SQL database for chemical probes: the compounds, the proteins they act
on, and the measurements linking the two.

```
database/schema.sql       the DDL, the single source of truth
database/probedb/         connect, insert, query
loader/                   read staging files, write them into the database
reference/                lookups that are the same for every source
staging/_template/        worked example of what a data source hands over
examples/                 populate_db.py, qc_db.py, explore_db.ipynb
tests.py
```

## 1. Clone

```bash
git clone git@github.com:SCHMM-MelMar/Summer-School-2026.git
cd Summer-School-2026
```

## 2. Install

```bash
# conda or mamba
mamba create -n db2 python=3.11 pandas loguru
mamba activate db2

# or a plain virtualenv
python -m venv .venv && source .venv/bin/activate
```

Install this repository editable:

```bash
pip install -e .
```

Check it worked:

```bash
python tests.py        # prints ok
```

## 3. Prepare your data

A data source hands over one directory of four flat files, `.tsv` or `.csv`.

```bash
cp -r staging/_template staging/my_source
```

| file                | columns                                                                                                                                                                                                                     |
| ------------------- | --------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------- |
| `compound.tsv`    | `inchikey` · `smiles` · `chembl_id`· `name`                                                                                                                                                                                 |
| `target.tsv`      | `target_key` · `type` · `name`                                                                                                                                                                                              |
| `uniprot.tsv`     | `uniprot_id` · `target_key` · `hgnc` · `species`                                                                                                                                                                         |
| `bioactivity.tsv` | `inchikey` · `target_key` · `moa` · `bioactivity_type` · `relation` · `value` · `unit` · `assay_type` · `assay_description` · `cell_line` · `concentration` · `concentration_unit` · `source_db` · `source` · `source_xref` · `xref_id` |

Only `inchikey`, `target_key`, `value` and `unit` have to be there in
`bioactivity.tsv`, and `inchikey` in `compound.tsv`. Any other column can be
left out of the header entirely and it loads as empty.

There is no fifth file for provenance: the directory is the set. Loading
`staging/reFRAME` records a compound set called `reFRAME` and puts every
compound in it, so `db.compounds(set="reFRAME")` works with nothing extra to
maintain. A compound arriving from three directories is one compound row and
three memberships.

The rules, in short:

- `inchikey` is the compound key everywhere.
- `target_key` is whatever string your source uses. It is a join key inside your
  directory and never reaches the database. Several `uniprot` rows sharing one
  `target_key` declare a complex or a family, and the same accession may appear
  under more than one key.
- `type` on a target is one of the four below. It is the only column with a
  fixed vocabulary, and the database enforces it. `protein`, `complex`, `family`, `ppi`
- `unit` says what the number means. Nothing is converted on load, so a source
  reporting pIC50 on -log(M) stays that way.
- `relation` should always be written.

The four provenance columns say where one number came from:

- `source_db` is the resource, `opnMe` or `Probes & Drugs` or `literature`. If
  you leave the column out it falls back to the directory name.
- `source` is the record inside it, a paper, an internal report, a release.
- `source_xref` identifies the single measurement in that record, so somebody
  can go back to the exact row.
- `xref_id` is the prefix that makes `source_xref` a link. It is a property of
  the source, so write the same value on every row of it. Whether there is one
  depends on the raw data: for a paper it is `https://doi.org/` with the DOI in
  `source_xref`, for an internal report there is nothing to resolve and you
  leave it empty.


Check a directory before handing it over:

```python
from loader import validate

validate("staging/my_source")     # [] means it is fine
```

The loader skips any directory whose name starts with `_`, so the template is
never loaded as data.

## 4. Populate

```python
from probedb import ProbeDB
from loader import load_all

db = ProbeDB("probe.db", create=True)   # a file on disk
load_all(db)                            # every directory under staging/
```

Pass `":memory:"` instead of a path to build the database in RAM. Nothing is
written to disk and it disappears when you are done, which is what `tests.py`
uses.

`create=True` builds the tables, so it fails on a database that already has
them. Delete the file first for a clean rebuild. To open one that already
exists, leave it out:

```python
db = ProbeDB("probe.db")     # read what is already there
```

Assembling all of `staging/` takes about two minutes, so it is worth keeping
the file. `explore_db.ipynb` opens it rather than rebuilding.

One source at a time:

```python
from loader import load

load(db, "staging/my_source", source="my_source")
```

`source` names the compound set, and stands in for any measurement whose row
does not say which `source_db` it came from. `load_all` passes the directory
name, so the set is called after the directory.

From the shell:

```bash
python examples/populate_db.py                 # refuse a directory that does not validate
python examples/populate_db.py --lenient       # load what each one can represent
```

`--lenient` is for the merged database while sources are still being cleaned
up. It loads what a directory can represent instead of refusing the whole
thing, and counts everything it dropped:

```
reFRAME  {'compounds': 8513, 'targets': 5609, 'sets': 1, 'bioactivities': 435266}
reFRAME  dropped {'targets_skipped': 1642, 'bioactivities_skipped': 1005973}
```

Loading the same directory twice does not duplicate anything. Compounds,
targets and sources are matched on their natural keys, and a measurement counts
as one already there when the compound, target, source, endpoint, operator,
value and unit all match. The report says how many rows were skipped:

```python
load(db, "staging/my_source")
# {'bioactivities': 14, 'duplicates_skipped': 0, ...}
```

Two papers reporting different numbers for the same pair stay two rows, which is
the point. Only exact repeats are collapsed.

```python
db.counts()
```

```
              table  rows
           compound     3
             chembl     3
       compound_set     1
compound_set_member     3
            uniprot     8
            target     6
    target_uniprot     9
bioactivity_source    10
 bioactivity_group     8
       bioactivity    13
```

### Checking a merged database

`populate_db.py` says what each directory contributed. `qc_db.py` asks the
opposite question: now that six sources are in one file, did they actually
merge, or is it six databases sharing a file?

```bash
python examples/qc_db.py probe.db
```

It splits what it finds in two. Checks are things that would mean the database
is wrong, and it exits non-zero if any fail:

```
    ok    one protein is one target however many sources report it
              0 accessions ended up on two rows
    ok    compounds are shared between sets
              2631 in more than one set
    FAIL  no log endpoint carries a molar unit
              46 rows
```

Notes are things worth knowing that belong to the source data rather than to
the merge, so they are printed and not failed:

```
    measurements with no unit: 6374
             source_db    n  no_type
        Probes & Drugs 6224     6224
```

## 5. Fetching data

Everything returns a pandas DataFrame. The outputs below are from the template.

### Search for a compound by name

```python
db.find("BI-")
```

```
                   inchikey    name                                     smiles
XQVVPGYIWAGRNI-JOCHJYFZSA-N BI-2536 N1([C@@H](C(=O)N(c2c1nc(nc2)Nc3c(cc(cc3)...
```

`db.bioactivities` and `db.compound_key` accept a name, an InChIKey or a ChEMBL
id, so all three of these are the same compound:

```python
db.compound_key("Olaparib")
db.compound_key("CHEMBL521686")
db.compound_key("FDLYAMZZIXQODN-UHFFFAOYSA-N")
```

### One row per compound, and which set it came from

`db.compounds()` is the compound index: one row per structure whatever it is a
member of, with the sets it belongs to spelled out.

```python
db.compounds()[["name", "n_sets", "sets", "n_targets"]]
```

```
    name  n_sets     sets  n_targets
 BI-2536       1 template          3
Olaparib       1 template          3
 (+)-JQ1       1 template          2
```

Filtering keeps that shape, so a compound in five libraries is still one row
and still shows all five:

```python
db.compounds(set="reFRAME")        # by set, which is the directory it came in from
db.compounds(category="library")   # by kind of collection
```

The other direction, from a compound to what claims it:

```python
db.compound_sets("BI-2536")
```

```
 set_id     name category source_db description
      1 template  library  template        None
```

And `db.sets()` is the catalogue, with how big each one is:

```python
db.sets()
```

```
 set_id     name category source_db  compounds
      1 template  library  template          3
```

Membership is a link table, not a column on the compound, so a compound is in
as many sets as claim it and the overlap between two libraries is a join:

```python
db.read("""
  SELECT a.name AS set_a, b.name AS set_b, COUNT(*) AS shared
    FROM compound_set_member ma JOIN compound_set a ON a.set_id = ma.set_id
    JOIN compound_set_member mb ON mb.inchikey = ma.inchikey
    JOIN compound_set b ON b.set_id = mb.set_id
   WHERE a.set_id < b.set_id GROUP BY a.set_id, b.set_id ORDER BY shared DESC
""")
```

If you want one row per membership instead, to merge a frame of your own onto,
that is the `compound_flat` view.

### Protein families

`uniprot.superfamily` is UniProt's own classification of a protein, filled in
from `reference/uniprot_protein_families.tsv` rather than from any one source,
so every source agrees on it. It is the whole string UniProt gives, outermost
group first:

```
P01116  KRAS   Small GTPase superfamily, Ras family
P53350  PLK1   Protein kinase superfamily, Ser/Thr protein kinase family, CDC5/Polo subfamily
```

`db.families()` splits that into levels and counts what sits under each, so you
can ask for either the superfamily or the family by name:

```python
db.families(like="kinase")
```

```
                            family  proteins  targets
        Protein kinase superfamily       249      474
         Tyr protein kinase family        58      120
 AGC Ser/Thr protein kinase family        37       61
```

Then `family=` on the three query methods:

```python
db.targets(family="Ras family")             # the targets in it, with all their members
db.compounds(family="Ras family")           # what has been measured against them
db.bioactivities(family="Ras family")       # the measurements themselves
```

Asking for `Protein kinase superfamily` gets everything under it; asking for
`Tyr protein kinase family` gets only that level.

Two things to know before trusting an answer.

**It is measured against, not active against.** Deciding what counts as active
needs a threshold and a unit, and the database picks neither, so a counter
screen at `Kd > 30 uM` comes back too. Filter on `value` and `relation`
yourself if you want hits.

**A family answer is only as good as the target rows underneath it.** A target
of type `protein` should have one accession. 258 of them have more, because a
source filed a PROTAC screen or an ortholog group as a single protein, and each
of those accessions drags the target into its own family. `qc_db.py` lists the
worst offenders. Check `db.targets(family=...)` before quoting a number from
`db.compounds(family=...)`.

Note that this is a property of a *protein*, and separate from a target of type
`family`, which is our own curated grouping like `PARP 1, 2 and 3`. A complex
has as many classifications as it has members.

### Everything measured for one compound

```python
df = db.bioactivities(compound="BI-2536")
df[["target", "moa", "bioactivity_type", "relation", "value", "unit",
    "assay_type", "cell_line", "source_db"]]
```

```
                              target       moa bioactivity_type relation    value    unit  assay_type cell_line      source_db
Serine/threonine-protein kinase PLK1 inhibitor             IC50        =     0.83      nM biochemical       NaN          opnMe
Serine/threonine-protein kinase PLK1 inhibitor            pIC50        =    10.08 -log(M) biochemical       NaN Probes & Drugs
Serine/threonine-protein kinase PLK1 inhibitor             IC50        =     1.10      nM biochemical       NaN       in-house
Serine/threonine-protein kinase PLK1 inhibitor             EC50        =    12.00      nM        cell      HeLa       in-house
    Bromodomain-containing protein 4 inhibitor             IC50        =     1.20      nM biochemical       NaN     literature
 Cyclin-dependent kinase 1/cyclin B1                       IC50        > 10000.00      nM biochemical       NaN     literature
```

Four numbers for PLK1 from three databases, one of them on a different scale,
plus a cell assay and a counter screen saying the compound does not hit
CDK1/cyclin B1. Nothing is converted, so the unit column is what each number
means.

### Everything measured against one target

Accepts a UniProt accession, an HGNC symbol or a target id.

```python
db.bioactivities(target="BRD4")[
    ["compound", "bioactivity_type", "relation", "value", "unit", "source"]]
```

```
compound bioactivity_type relation  value unit        source
 BI-2536             IC50        =    1.2   nM PMID:32088495
 (+)-JQ1             EC50        =  144.5   nM PMID:31303996
```

### One accession, several targets

An accession can belong to a protein and to the groups containing it, so this
returns a list.

```python
db.targets_for("P09874")          # -> [4, 5]
```

```
4  protein  Poly [ADP-ribose] polymerase 1
5  family   PARP 1, 2 and 3
```

Those are two different targets and two different sets of measurements.

### Targets together with their accessions

`db.table("target")` and `db.table("uniprot")` have no column in common, because
the link between them is the third table `target_uniprot`. `db.targets()` does
both joins and hands back one flat frame, one row per target and accession pair:

```python
t = db.targets()
t[t.type != "protein"]
```

```
 target_id    type                                name uniprot_id  hgnc      species entrez_gene
         5  family                     PARP 1, 2 and 3     P09874 PARP1 Homo sapiens        None
         5  family                     PARP 1, 2 and 3     Q9UGN5 PARP2 Homo sapiens        None
         5  family                     PARP 1, 2 and 3     Q9Y6F1 PARP3 Homo sapiens        None
         6 complex Cyclin-dependent kinase 1/cyclin B1     P06493  CDK1 Homo sapiens        None
         6 complex Cyclin-dependent kinase 1/cyclin B1     P14635 CCNB1 Homo sapiens        None
```

It takes the same argument as `bioactivities(target=...)`, an accession, an HGNC
symbol or a target id, so `db.targets("PARP1")` gives the protein and the family
containing it. The same frame is the `target_flat` view if you are writing SQL,
and `target_id` joins it straight onto `db.bioactivities()`.

If you would rather keep the raw tables, the merge is two steps:

```python
(db.table("target")
   .merge(db.table("target_uniprot"), on="target_id")
   .merge(db.table("uniprot"), on="uniprot_id"))
```

To find a target by its members instead:

```python
db.find_target("family", ["Q9Y6F1", "P09874", "Q9UGN5"])   # order does not matter
```

### Where a number came from

`db.sources()` lists every source with how much it contributed.

```python
db.sources()
```

```
 source_id      source_db               source          xref_id  measurements
         3       in-house assay report 2024-03              NaN             3
         7     literature        PMID:35007061 https://doi.org/             2
         2 Probes & Drugs        PMID:33539089              NaN             1
         5     literature        PMID:17291758 https://doi.org/             1
        10     literature        PMID:23473053 https://doi.org/             1
         6     literature        PMID:31303996 https://doi.org/             1
```

`source_id` is on every measurement, so this joins onto `db.bioactivities()`
directly. Where the source has an `xref_id`, `bioactivities()` builds the link
for you in `source_url`:

```python
act = db.bioactivities()
act[act.source_url.notna()][
    ["compound", "target", "value", "unit", "source_db", "source", "source_url"]]
```

```
compound                              target   value unit  source_db        source                                   source_url
 BI-2536    Bromodomain-containing protein 4     1.2   nM literature PMID:32088495 https://doi.org/10.1016/j.ejmech.2020.112152
 BI-2536 Cyclin-dependent kinase 1/cyclin B1 10000.0   nM literature PMID:17291758    https://doi.org/10.1016/j.cub.2006.12.037
 (+)-JQ1    Bromodomain-containing protein 4   144.5   nM literature PMID:31303996           https://doi.org/10.1039/c8md00412a
 (+)-JQ1    Bromodomain-containing protein 2   120.2   nM literature PMID:35007061 https://doi.org/10.1021/acs.jmedchem.1c01779
```

`source_url` is empty when the source has nothing resolvable, which is the case
for the in-house rows and for a Probes & Drugs activity id. The identifier is
still in `source_xref` either way.

`bioactivity` itself stores `source_id`, an integer, because the name of a
source belongs in one place and not on every one of its measurements. If you
want the names without writing the join, read `bioactivity_flat` instead of
`bioactivity`. It is the same rows with the compound, the target and the source
spelled out, and it is what `db.bioactivities()` returns:

```python
db.read("SELECT source_db, COUNT(*) n FROM bioactivity_flat GROUP BY 1")
```

```
     source_db    n
Probes & Drugs 1684
         opnMe   75
```

### Assay context

How a number was measured sits on the measurement, so a biochemical IC50 and a
cell based EC50 for the same pair stay distinguishable.

```python
act = db.bioactivities()
act[act.assay_type.notna()][
    ["compound", "target", "bioactivity_type", "value", "unit", "assay_type",
     "assay_description", "cell_line"]].head()
```

```
compound                               target bioactivity_type  value    unit  assay_type     assay_description cell_line
 BI-2536 Serine/threonine-protein kinase PLK1             IC50   0.83      nM biochemical kinase activity assay       NaN
 BI-2536 Serine/threonine-protein kinase PLK1            pIC50  10.08 -log(M) biochemical                   NaN       NaN
 BI-2536 Serine/threonine-protein kinase PLK1             IC50   1.10      nM biochemical kinase activity assay       NaN
```


### Fixed vocabularies

The closed vocabularies live in `CHECK` constraints in `database/schema.sql`
and nowhere else. The loader reads them back out of that file, so there is one
definition and it cannot drift:

```python
from probedb.schema import vocabulary

vocabulary("type")       # {'protein', 'complex', 'family', 'ppi'}
vocabulary("relation")   # {'=', '>', '<', '>=', '<=', '~'}
vocabulary("category")   # what kind of collection a compound set is
```

Adding a value means editing the `CHECK` in `schema.sql`. Everything else
follows.

## 6. The schema

| table                  | holds                                                        |
| ---------------------- | ------------------------------------------------------------ |
| `compound`           | one row per structure, keyed on InChIKey, searchable by name |
| `chembl`             | ChEMBL id of a compound, keyed on the id                     |
| `compound_set`       | a collection, one per staging directory                      |
| `compound_set_member`| which compounds are in it                                    |
| `uniprot`            | accession, HGNC symbol, species, UniProt family classification |
| `target`             | one row per target, with a type                              |
| `target_uniprot`     | which accessions make up a target                            |
| `bioactivity_source` | where the numbers came from, and how to resolve an xref      |
| `bioactivity_group`  | compound, target and MoA, the unit you would aggregate over  |
| `bioactivity`        | one row per reported measurement, with how it was measured   |
| `target_flat`        | view, one row per target and accession pair                  |
| `compound_flat`      | view, one row per compound and set                           |
| `bioactivity_flat`   | view, every measurement with its compound, target and source named |

A few things worth knowing before you query.

**The InChIKey is the compound key.** There is no surrogate compound id: every
source is matched on the InChIKey anyway, so a second identifier would only be
another name for the same thing. `bioactivity` and `bioactivity_group` carry it
directly, and `chembl` maps it to ChEMBL ids.

**A target is not an accession.** `target` says what kind of thing it is,
`target_uniprot` says which proteins make it up. Identity is the member set, so
the same group reported by two sources in two orders lands on one row, while
PARP1 and PARP 1, 2 and 3 stay distinct. A single protein is a target with one
member, so nothing has to special case it.

**Units are kept as reported.** Nothing is converted anywhere, and nothing is
filtered on plausibility. `relation` carries the operator, so a censored
`IC50 > 10 uM` is stored without being mistakable for a potency, and any
comparison has to name the unit it is comparing on.

**Measurements are never merged.** The same compound, target and assay measured
three times stays three rows. Averaging on load would hide the spread and lose
the per source attribution.

**A set is a directory, and membership is a table.** The same compound arriving
from five sources is one `compound` row and five `compound_set_member` rows, so
where it came from survives the merge and the overlap between two sources is a
join. A column on `compound` could not hold five answers.
