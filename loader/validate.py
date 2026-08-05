import csv
import re
from pathlib import Path

from probedb.db import INCHIKEY
from probedb.schema import vocabulary

FILES = {
    "compound": (["inchikey"], ["smiles", "chembl_id", "name"]),
    "target": (["target_key"], ["type", "name"]),
    "uniprot": (["uniprot_id", "target_key"], ["hgnc", "species", "entrez_gene"]),
    "bioactivity": (
        ["inchikey", "target_key", "value", "unit"],
        [
            "moa",
            "bioactivity_type",
            "relation",
            "source_db",
            "source",
            "source_xref",
            "assay_type",
            "assay_description",
            "cell_line",
            "concentration",
            "concentration_unit",
        ],
    ),
}

# the vocabularies are defined once, in the CHECK constraints in schema.sql
RELATIONS = vocabulary("relation")
TARGET_TYPES = vocabulary("type")
CHEMBL = re.compile(r"^CHEMBL[0-9]+$", re.I)


def read(directory, name):
    for suffix, delimiter in ((".tsv", "\t"), (".csv", ",")):
        path = Path(directory) / (name + suffix)
        if path.exists():
            with open(path, newline="") as handle:
                return [
                    {k.strip(): (v or "").strip() for k, v in row.items()}
                    for row in csv.DictReader(handle, delimiter=delimiter)
                ]
    return []


def validate(directory):
    problems = []
    tables = {name: read(directory, name) for name in FILES}

    for name, (required, optional) in FILES.items():
        rows = tables[name]
        if not rows:
            if name in ("compound", "bioactivity"):
                problems.append(f"{name}: file missing or empty")
            continue
        missing = set(required) - set(rows[0])
        if missing:
            problems.append(f"{name}: missing column(s) {sorted(missing)}")
        unknown = set(rows[0]) - set(required) - set(optional)
        if unknown:
            problems.append(f"{name}: unknown column(s) {sorted(unknown)}, ignored")

    # the InChIKey is the compound primary key, and a ChEMBL id points at one
    # compound, so both have to be consistent inside the directory
    keys, owner = set(), {}
    for line, row in enumerate(tables["compound"], 2):
        inchikey = row.get("inchikey", "").upper()
        if not INCHIKEY.match(inchikey):
            problems.append(
                f"compound line {line}: {row.get('inchikey')!r} "
                f"is not a standard InChIKey"
            )
            continue
        if inchikey in keys:
            problems.append(
                f"compound line {line}: {inchikey} appears more than once, "
                f"only the first row is used, kept anyway"
            )
        keys.add(inchikey)

        for chembl_id in (row.get("chembl_id") or "").replace(";", "|").split("|"):
            chembl_id = chembl_id.strip().upper()
            if not chembl_id:
                continue
            if not CHEMBL.match(chembl_id):
                problems.append(
                    f"compound line {line}: {chembl_id!r} does not look like a "
                    f"ChEMBL id, kept anyway"
                )
            if owner.setdefault(chembl_id, inchikey) != inchikey:
                problems.append(
                    f"compound line {line}: {chembl_id} is already used by "
                    f"{owner[chembl_id]}, a ChEMBL id belongs to one compound"
                )

    declared = {}
    for line, row in enumerate(tables["target"], 2):
        if row.get("type") and row["type"] not in TARGET_TYPES:
            problems.append(
                f"target line {line}: type {row['type']!r} is not one of "
                f"{sorted(TARGET_TYPES)}"
            )
        key = row.get("target_key")
        if key in declared and declared[key] != (row.get("type"), row.get("name")):
            problems.append(
                f"target line {line}: {key!r} is declared twice with a different "
                f"type or name"
            )
        declared[key] = (row.get("type"), row.get("name"))

    target_keys = {r.get("target_key") for r in tables["target"]}
    for line, row in enumerate(tables["uniprot"], 2):
        if row.get("target_key") not in target_keys:
            problems.append(
                f"uniprot line {line}: unknown target_key " f"{row.get('target_key')!r}"
            )

    for line, row in enumerate(tables["bioactivity"], 2):
        if row.get("inchikey", "").upper() not in keys:
            problems.append(
                f"bioactivity line {line}: unknown inchikey " f"{row.get('inchikey')!r}"
            )
        if row.get("target_key") not in target_keys:
            problems.append(
                f"bioactivity line {line}: unknown target_key "
                f"{row.get('target_key')!r}"
            )
        if not row.get("unit"):
            problems.append(f"bioactivity line {line}: no unit, kept anyway")
        if not row.get("relation"):
            problems.append(f"bioactivity line {line}: no relation, kept anyway")
        if row.get("relation") and row["relation"] not in RELATIONS:
            problems.append(
                f"bioactivity line {line}: bad relation " f"{row['relation']!r}"
            )
        for column in ("value", "concentration"):
            try:
                if row.get(column):
                    float(row[column])
            except ValueError:
                problems.append(
                    f"bioactivity line {line}: {column} "
                    f"{row[column]!r} is not numeric"
                )

    return problems
