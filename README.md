# DB

A small SQL database for chemical probes: the compounds, the proteins they act
on, and the measurements linking the two.

```
database/schema.sql       the DDL, the single source of truth
database/probedb/         connect, insert, query
loader/                   read staging files, write them into the database
staging/_template/        worked example of what a data source hands over
chemicalprobes.org/       the same thing done from a real source, start at its README
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
| `bioactivity.tsv` | `inchikey` · `target_key` · `moa` · `bioactivity_type` · `relation` · `value` · `unit` · `assay_type` · `assay_description` · `cell_line` · `concentration` · `concentration_unit` · `source_db` · `source` · `source_xref` · `xref_id` |

Only `inchikey`, `target_key`, `value` and `unit` have to be there in
`bioactivity.tsv`, and `inchikey` in `compound.tsv`. Any other column can be
left out of the header entirely and it loads as empty.

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
bioactivity_source    10
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
| `chembl`             | ChEMBL id of a compound, keyed on the id                     |
| `uniprot`            | accession, HGNC symbol, species                              |
| `target`             | one row per target, with a type                              |
| `target_uniprot`     | which accessions make up a target                            |
| `bioactivity_source` | where the numbers came from, and how to resolve an xref      |
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
