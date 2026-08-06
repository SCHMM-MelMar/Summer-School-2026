import re
import sqlite3

import pandas as pd

from . import schema

INCHIKEY = re.compile(r"^[A-Z]{14}-[A-Z]{10}-[A-Z]$")


class ProbeDB:

    def __init__(self, path=":memory:", create=False):
        self.conn = sqlite3.connect(str(path))
        self.conn.execute("PRAGMA foreign_keys = ON")
        if create:
            self.create()

    def create(self):
        for statement in schema.statements():
            self.conn.execute(statement)
        self.conn.commit()

    def commit(self):
        self.conn.commit()

    def close(self):
        self.conn.close()

    # reading

    def read(self, sql, *params):
        return pd.read_sql_query(sql, self.conn, params=params)

    def one(self, sql, *params):
        row = self.conn.execute(sql, params).fetchone()
        return None if row is None else row[0]

    def table(self, name):
        return self.read(f"SELECT * FROM {name}")

    def counts(self):
        rows = [
            (name, self.one(f"SELECT COUNT(*) FROM {name}")) for name in schema.TABLES
        ]
        return pd.DataFrame(rows, columns=["table", "rows"])

    def targets(self, key=None):
        # target and uniprot share no column, target_uniprot is the link.
        # one row per target and accession, so a complex comes back as
        # several rows with the same target_id
        sql = "SELECT * FROM target_flat"
        params = []
        if key is not None:
            ids = self.targets_for(key)
            sql += f" WHERE target_id IN ({','.join('?' * len(ids))})"
            params = ids
        return self.read(sql + " ORDER BY target_id, uniprot_id", *params)

    def sources(self):
        return self.read(
            """
            SELECT s.source_id, s.source_db, s.source, s.xref_id,
                   COUNT(b.id) AS measurements
              FROM bioactivity_source s
              LEFT JOIN bioactivity b ON b.source_id = s.source_id
             GROUP BY s.source_id
             ORDER BY measurements DESC, s.source_db, s.source
        """
        )

    def bioactivities(self, compound=None, target=None):
        where, params = [], []
        if compound is not None:
            where.append("b.inchikey = ?")
            params.append(self.compound_key(compound))
        if target is not None:
            ids = self.targets_for(target)
            where.append(f"b.target_id IN ({','.join('?' * len(ids))})")
            params += ids
        sql = """
            SELECT b.id, b.inchikey, c.name AS compound, b.target_id,
                   t.type AS target_type, t.name AS target, b.moa,
                   b.bioactivity_type, b.relation, b.value, b.unit,
                   b.assay_type, b.assay_description, b.cell_line,
                   b.concentration, b.concentration_unit,
                   b.source_id, s.source_db, s.source, b.source_xref,
                   CASE WHEN s.xref_id IS NOT NULL AND b.source_xref IS NOT NULL
                        THEN s.xref_id || b.source_xref END AS source_url
              FROM bioactivity b
              JOIN compound c ON c.inchikey = b.inchikey
              JOIN target t ON t.target_id = b.target_id
              LEFT JOIN bioactivity_source s ON s.source_id = b.source_id
        """
        if where:
            sql += " WHERE " + " AND ".join(where)
        return self.read(sql + " ORDER BY b.id", *params)

    # lookup

    def compound_key(self, key):
        key = str(key).strip()
        found = self.one(
            "SELECT inchikey FROM compound WHERE inchikey = ?", key.upper()
        )
        if found is not None:
            return found
        found = self.one("SELECT inchikey FROM chembl WHERE chembl_id = ?", key.upper())
        if found is not None:
            return found
        return self.one(
            "SELECT inchikey FROM compound WHERE name = ? COLLATE NOCASE", key
        )

    def find(self, name):
        return self.read(
            "SELECT inchikey, name, smiles FROM compound "
            "WHERE name LIKE ? ORDER BY name",
            f"%{name.strip()}%",
        )

    def targets_for(self, key):
        # an accession belongs to several targets: P09874 is PARP1 on its own
        # and a member of the PARP1/2/3 complex
        if isinstance(key, int):
            return [key]
        key = str(key).strip().upper()
        rows = self.conn.execute(
            "SELECT DISTINCT tu.target_id FROM target_uniprot tu "
            "  JOIN uniprot u ON u.uniprot_id = tu.uniprot_id "
            " WHERE tu.uniprot_id = ? OR u.hgnc = ? ORDER BY tu.target_id",
            (key, key),
        ).fetchall()
        return [r[0] for r in rows]

    def find_target(self, type, accessions):
        accessions = sorted({a.strip().upper() for a in accessions if a})
        if not accessions:
            return None
        marks = ",".join("?" * len(accessions))
        return self.one(
            f"""
            SELECT tu.target_id FROM target_uniprot tu
              JOIN target t ON t.target_id = tu.target_id
             WHERE t.type = ?
             GROUP BY tu.target_id
            HAVING COUNT(*) = ?
               AND SUM(CASE WHEN tu.uniprot_id IN ({marks}) THEN 1 ELSE 0 END) = ?
        """,
            type,
            len(accessions),
            *accessions,
            len(accessions),
        )

    # writing

    def insert(self, table, **values):
        columns = ", ".join(values)
        marks = ", ".join("?" * len(values))
        cursor = self.conn.execute(
            f"INSERT INTO {table} ({columns}) VALUES ({marks})", tuple(values.values())
        )
        return cursor.lastrowid

    def add_compound(self, inchikey, smiles=None, name=None):
        inchikey = (inchikey or "").strip().upper()
        if not INCHIKEY.match(inchikey):
            raise ValueError(f"not a standard InChIKey: {inchikey!r}")
        row = self.conn.execute(
            "SELECT smiles, name FROM compound WHERE inchikey = ?", (inchikey,)
        ).fetchone()
        if row is None:
            self.insert(
                "compound", inchikey=inchikey, smiles=smiles or None, name=name or None
            )
            return inchikey
        if smiles and not row[0]:
            self.conn.execute(
                "UPDATE compound SET smiles = ? WHERE inchikey = ?", (smiles, inchikey)
            )
        if name and not row[1]:
            self.conn.execute(
                "UPDATE compound SET name = ? WHERE inchikey = ?", (name, inchikey)
            )
        return inchikey

    def add_chembl(self, inchikey, chembl_id):
        inchikey, chembl_id = inchikey.strip().upper(), chembl_id.strip().upper()
        if (
            self.one("SELECT chembl_id FROM chembl WHERE chembl_id = ?", chembl_id)
            is None
        ):
            self.insert("chembl", chembl_id=chembl_id, inchikey=inchikey)
        return chembl_id

    def add_uniprot(self, uniprot_id, hgnc=None, species=None, entrez_gene=None):
        uniprot_id = uniprot_id.strip().upper()
        if (
            self.one("SELECT uniprot_id FROM uniprot WHERE uniprot_id = ?", uniprot_id)
            is None
        ):
            self.insert(
                "uniprot",
                uniprot_id=uniprot_id,
                hgnc=hgnc or None,
                species=species or None,
                entrez_gene=entrez_gene or None,
            )
        return uniprot_id

    def add_target(self, type, name=None, uniprots=()):
        uniprots = [u if isinstance(u, (list, tuple)) else (u,) for u in uniprots]
        accessions = [u[0].strip().upper() for u in uniprots if u[0]]

        if accessions:
            found = self.find_target(type, accessions)
        else:
            found = self.one(
                "SELECT target_id FROM target WHERE type = ? AND name = ?", type, name
            )
        if found is not None:
            return found

        target_id = self.insert("target", type=type, name=name or None)
        for uniprot in uniprots:
            self.add_uniprot(*uniprot[:4])
        for accession in sorted(set(accessions)):
            self.insert("target_uniprot", target_id=target_id, uniprot_id=accession)
        return target_id

    def add_source(self, source_db, source=None, xref_id=None):
        found = self.one(
            "SELECT source_id FROM bioactivity_source "
            "WHERE source_db = ? AND source IS ?",
            source_db,
            source,
        )
        if found is None:
            return self.insert(
                "bioactivity_source",
                source_db=source_db,
                source=source or None,
                xref_id=xref_id or None,
            )
        if xref_id:
            self.conn.execute(
                "UPDATE bioactivity_source SET xref_id = ? "
                "WHERE source_id = ? AND xref_id IS NULL",
                (xref_id, found),
            )
        return found

    def add_group(self, inchikey, target_id, moa=""):
        found = self.one(
            "SELECT id FROM bioactivity_group WHERE inchikey = ? "
            "AND target_id = ? AND moa = ?",
            inchikey,
            target_id,
            moa,
        )
        return (
            found
            if found is not None
            else self.insert(
                "bioactivity_group", inchikey=inchikey, target_id=target_id, moa=moa
            )
        )

    def add_bioactivity(
        self,
        inchikey,
        target_id,
        moa="",
        bioactivity_type=None,
        relation=None,
        value=None,
        unit=None,
        assay_type=None,
        assay_description=None,
        cell_line=None,
        concentration=None,
        concentration_unit=None,
        source_id=None,
        source_xref=None,
    ):
        moa = moa or ""
        self.add_group(inchikey, target_id, moa)
        return self.insert(
            "bioactivity",
            inchikey=inchikey,
            target_id=target_id,
            moa=moa,
            bioactivity_type=bioactivity_type,
            relation=relation or None,
            value=None if value in (None, "") else float(value),
            unit=unit,
            assay_type=assay_type or None,
            assay_description=assay_description or None,
            cell_line=cell_line or None,
            concentration=None if concentration in (None, "") else float(concentration),
            concentration_unit=concentration_unit or None,
            source_id=source_id,
            source_xref=source_xref or None,
        )
