import csv
import io
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
            "xref_id",
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

# A target is identified by the set of accessions it is made of, so a
# placeholder in this column is not a missing value, it is a wrong one: every
# target carrying it looks like the same target and they merge on load.
UNIPROT = re.compile(
    r"^([OPQ][0-9][A-Z0-9]{3}[0-9]|[A-NR-Z][0-9]([A-Z][A-Z0-9]{2}[0-9]){1,2})$"
)


def text(path):
    """The file's contents, and the encoding it turned out to be written in."""
    # a spreadsheet saved out of Excel is often cp1252, and one micro sign in a
    # free text column is not a reason to refuse a whole directory
    raw = path.read_bytes()
    for encoding in ("utf-8", "cp1252"):
        try:
            return raw.decode(encoding), encoding
        except UnicodeDecodeError:
            pass
    return raw.decode("utf-8", errors="replace"), "utf-8 with bytes replaced"


DELIMITERS = {"tab": "\t", "comma": ",", "semicolon": ";"}


def delimiter_of(body, expected):
    """What actually separates the header, whatever the file is named."""
    header = body.split("\n", 1)[0]
    if header.count(DELIMITERS[expected]):
        return expected
    found = max(DELIMITERS, key=lambda name: header.count(DELIMITERS[name]))
    return found if header.count(DELIMITERS[found]) else expected


def read(directory, name, layout=None):
    for suffix, expected in ((".tsv", "tab"), (".csv", "comma")):
        path = Path(directory) / (name + suffix)
        if path.exists():
            body, encoding = text(path)
            delimiter = delimiter_of(body, expected)
            if layout is not None:
                layout[path.name] = (encoding, expected, delimiter)
            return [
                {k.strip(): (v or "").strip() for k, v in row.items()}
                for row in csv.DictReader(
                    io.StringIO(body, newline=""), delimiter=DELIMITERS[delimiter]
                )
            ]
    return []


def validate(directory):
    problems = []
    layout = {}
    tables = {name: read(directory, name, layout) for name in FILES}

    for filename, (encoding, expected, delimiter) in sorted(layout.items()):
        if encoding != "utf-8":
            problems.append(
                f"{filename}: not UTF-8, read as {encoding}. Re-save it as "
                f"UTF-8 or the characters that are not plain ASCII will be "
                f"wrong, kept anyway"
            )
        if delimiter != expected:
            problems.append(
                f"{filename}: named for {expected} separated but is "
                f"{delimiter} separated, read as {delimiter}. Rename the file "
                f"so the two agree, kept anyway"
            )

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
            # the InChIKey is the join key across all four files, so it cannot
            # be filled in later and it cannot be stood in for
            fix = (
                "derive it from the SMILES in your extractor"
                if row.get("smiles")
                else "a compound with no structure has no key and cannot be loaded"
            )
            problems.append(
                f"compound line {line}: {row.get('inchikey')!r} is not a "
                f"standard InChIKey, {fix}. Never substitute a placeholder, "
                f"every compound sharing one would merge into a single row. "
                f"skipped"
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
                f"uniprot line {line}: unknown target_key "
                f"{row.get('target_key')!r}, skipped"
            )
        accession = row.get("uniprot_id", "").upper()
        if not UNIPROT.match(accession):
            problems.append(
                f"uniprot line {line}: {row.get('uniprot_id')!r} is not a "
                f"UniProt accession, so the target keeps its name and loses "
                f"this member. Leave the row out rather than filling it in, "
                f"every target sharing a placeholder would merge into one. "
                f"skipped"
            )

    for line, row in enumerate(tables["bioactivity"], 2):
        if row.get("inchikey", "").upper() not in keys:
            problems.append(
                f"bioactivity line {line}: unknown inchikey "
                f"{row.get('inchikey')!r}, skipped"
            )
        if row.get("target_key") not in target_keys:
            problems.append(
                f"bioactivity line {line}: unknown target_key "
                f"{row.get('target_key')!r}, skipped"
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
