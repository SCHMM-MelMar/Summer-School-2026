# DB

A small SQL database for chemical probes: the compounds, the proteins they act
on, and the measurements linking the two.

```
database/schema.sql       the DDL, the single source of truth
database/probedb/         connect, insert, query
loader/                   read staging files, write them into the database
staging/_template/        worked example of what a data source hands over
examples/                 populate_db.py, explore_db.ipynb
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
| `bioactivity.tsv` | `inchikey` · `target_key` · `moa` · `bioactivity_type` · `relation` · `value` · `unit` · `assay_type` · `assay_description` · `cell_line` · `concentration` · `concentration_unit` · `source_db` · `source` · `source_xref` |

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
written to disk and it disappears when you are done, which is what the notebook
and `tests.py` use. Rebuilding from the staging files takes under a second, so
there is rarely a reason to keep a file around while exploring.

`create=True` builds the tables, so it fails on a database that already has
them. Delete the file first for a clean rebuild.

One source at a time:

```python
from loader import load

load(db, "staging/my_source", source="my_source")
```

From the shell:

```bash
python examples/populate_db.py
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
           uniprot     8
            target     6
    target_uniprot     9
bioactivity_source     9
 bioactivity_group     8
       bioactivity    13
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

### Everything measured for one compound

```python
df = db.bioactivities(compound="BI-2536")
df[["target", "moa", "bioactivity_type", "relation", "value", "unit",
    "assay_type", "cell_line", "source_db"]]
```

```
                              target       moa bioactivity_type relation     value    unit  assay_type cell_line      source_db
Serine/threonine-protein kinase PLK1 inhibitor             IC50        =     0.830      nM biochemical                   opnMe
Serine/threonine-protein kinase PLK1 inhibitor             IC50        =     1.100      nM biochemical                in-house
Serine/threonine-protein kinase PLK1 inhibitor             EC50        =    12.000      nM        cell      HeLa    in-house
Serine/threonine-protein kinase PLK1 inhibitor            pIC50        =     9.080 -log(M) biochemical          Probes & Drugs
    Bromodomain-containing protein 4 inhibitor             IC50        =     1.202      nM biochemical                  ChEMBL
 Cyclin-dependent kinase 1/cyclin B1                       IC50        > 10000.000      nM biochemical                  ChEMBL
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
compound bioactivity_type relation   value unit        source
 BI-2536             IC50        =   1.202   nM PMID:32088495
 (+)-JQ1             EC50        = 144.500   nM PMID:31303996
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

### What a group is made of

`target_flat` is a view with one row per target and accession pair.

```python
db.read("SELECT * FROM target_flat WHERE type <> 'protein' ORDER BY target_id")
```

```
 target_id    type                                name uniprot_id  hgnc      species
         5  family                     PARP 1, 2 and 3     P09874 PARP1 Homo sapiens
         5  family                     PARP 1, 2 and 3     Q9UGN5 PARP2 Homo sapiens
         5  family                     PARP 1, 2 and 3     Q9Y6F1 PARP3 Homo sapiens
         6 complex Cyclin-dependent kinase 1/cyclin B1     P06493  CDK1 Homo sapiens
         6 complex Cyclin-dependent kinase 1/cyclin B1     P14635 CCNB1 Homo sapiens
```

To find one by its members instead:

```python
db.find_target("family", ["Q9Y6F1", "P09874", "Q9UGN5"])   # order does not matter
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
compound                               target bioactivity_type  value unit  assay_type     assay_description cell_line
 BI-2536 Serine/threonine-protein kinase PLK1             IC50  0.830   nM biochemical kinase activity assay
 BI-2536 Serine/threonine-protein kinase PLK1             IC50  1.100   nM biochemical kinase activity assay
 BI-2536 Serine/threonine-protein kinase PLK1             EC50 12.000   nM        cell     antiproliferation      HeLa
```


### Fixed vocabularies

The two closed vocabularies live in `CHECK` constraints in `database/schema.sql`
and nowhere else. The loader reads them back out of that file, so there is one
definition and it cannot drift:

```python
from probedb.schema import vocabulary

vocabulary("type")       # {'protein', 'complex', 'family', 'ppi'}
vocabulary("relation")   # {'=', '>', '<', '>=', '<=', '~'}
```

Adding a value means editing the `CHECK` in `schema.sql`. Everything else
follows.

## 6. The schema

| table                  | holds                                                        |
| ---------------------- | ------------------------------------------------------------ |
| `compound`           | one row per structure, keyed on InChIKey, searchable by name |
| `chembl`             | InChIKey and SMILES of compounds that are in ChEMBL          |
| `uniprot`            | accession, HGNC symbol, species                              |
| `target`             | one row per target, with a type                              |
| `target_uniprot`     | which accessions make up a target                            |
| `bioactivity_source` | `source_db` and `source` pairs                           |
| `bioactivity_group`  | compound, target and MoA, the unit you would aggregate over  |
| `bioactivity`        | one row per reported measurement, with how it was measured   |
| `target_flat`        | view, one row per target and accession pair                  |

Three things worth knowing before you query.

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
