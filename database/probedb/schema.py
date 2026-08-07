import re
from pathlib import Path


SCHEMA_SQL = Path(__file__).resolve().parent.parent / "schema.sql"

TABLES = [
    "compound",
    "chembl",
    "compound_set",
    "compound_set_member",
    "uniprot",
    "target",
    "target_uniprot",
    "bioactivity_source",
    "bioactivity_group",
    "bioactivity",
]

VIEWS = ["bioactivity_flat", "compound_flat", "target_flat"]


def vocabulary(column):
    """The allowed values for a column, read from its CHECK constraint."""
    # \b so asking for `type` does not find `relation_type`
    pattern = rf"\b{column}\s+IN\s*\(([^)]*)\)"
    match = re.search(pattern, SCHEMA_SQL.read_text(), re.I)
    if not match:
        raise KeyError(f"no CHECK vocabulary for {column!r}")
    return {v.strip().strip("'\"") for v in match.group(1).split(",")}


def statements():
    if not SCHEMA_SQL.exists():
        raise FileNotFoundError(
            f"{SCHEMA_SQL} not found -- install with `pip install -e .`"
        )
    sql = re.sub(r"--[^\n]*", "", SCHEMA_SQL.read_text())
    sql = re.sub(
        r"\bSERIAL\s+PRIMARY\s+KEY\b",
        "INTEGER PRIMARY KEY AUTOINCREMENT",
        sql,
        flags=re.I,
    )
    return [s.strip() for s in sql.split(";") if s.strip()]
