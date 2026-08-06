"""Load this handoff bundle into a fresh SQLite database.

    python new_load.py probe.db

Reads new_schema.sql and every new_<table>.tsv beside it. Nothing else is
needed: no repository, no dependencies, no network.

sqlite3's own `.import` cannot be used for this. It writes an empty cell as the
empty string rather than NULL, so the first bioactivity row fails
`CHECK (relation IS NULL OR relation IN (...))` -- '' is neither. This turns a
blank into NULL, which is what the schema means by it.
"""

import csv
import re
import sqlite3
import sys
from pathlib import Path

HERE = Path(__file__).resolve().parent


def tables(schema_sql):
    """The tables in the order the file declares them, which is an order that
    satisfies the foreign keys: compound before chembl, target before
    target_uniprot."""
    return re.findall(r"CREATE TABLE (\w+)", schema_sql)


def load(db_path):
    schema_sql = (HERE / "new_schema.sql").read_text()
    connection = sqlite3.connect(db_path)
    connection.executescript(schema_sql.replace("SERIAL PRIMARY KEY",
                                                "INTEGER PRIMARY KEY AUTOINCREMENT"))
    connection.execute("PRAGMA foreign_keys = ON")

    loaded = {}
    for table in tables(schema_sql):
        path = HERE / f"new_{table}.tsv"
        if not path.exists():
            continue
        with open(path, newline="") as handle:
            rows = list(csv.DictReader(handle, delimiter="\t"))
        if not rows:
            loaded[table] = 0
            continue
        columns = list(rows[0])
        marks = ", ".join("?" * len(columns))
        connection.executemany(
            f"INSERT INTO {table} ({', '.join(columns)}) VALUES ({marks})",
            # a blank cell is NULL, not the empty string
            [[row[column] if row[column] != "" else None for column in columns]
             for row in rows],
        )
        loaded[table] = len(rows)
    connection.commit()

    broken = connection.execute("PRAGMA foreign_key_check").fetchall()
    if broken:
        raise SystemExit(f"foreign keys do not resolve: {broken[:5]}")
    return connection, loaded


if __name__ == "__main__":
    target = sys.argv[1] if len(sys.argv) > 1 else "probe.db"
    if Path(target).exists():
        raise SystemExit(f"{target} already exists, delete it or name another file")
    connection, loaded = load(target)
    width = max(len(name) for name in loaded)
    for name, count in loaded.items():
        print(f"{count:>7}  {name}")
    print(f"{sum(loaded.values()):>7}  rows in {target}, foreign keys resolve")
