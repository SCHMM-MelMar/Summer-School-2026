import re
import sqlite3

import pandas as pd

from . import schema

INCHIKEY = re.compile(r"^[A-Z]{14}-[A-Z]{10}-[A-Z]$")

# uniprot.superfamily is the whole classification, outermost group first:
# "Small GTPase superfamily, Ras family". this matches one level of it by name,
# so asking for the superfamily gets everything under it and asking for the
# family gets just that family. LIKE is case insensitive in SQLite.
FAMILY_MATCH = "(', ' || u.superfamily || ', ') LIKE ('%, ' || ? || ', %')"


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

    def compounds(self, set=None, category=None, family=None):
        # one row per compound, whatever it is a member of. the sets are rolled
        # up into a readable column here and stay a link table underneath, so
        # filtering by one set still shows you the other sets it belongs to
        member, params = [], []
        if set is not None:
            member.append("s.name = ?")
            params.append(set)
        if category is not None:
            member.append("s.category = ?")
            params.append(category)

        where = []
        if member:
            where.append(
                f"""c.inchikey IN (
                    SELECT m.inchikey FROM compound_set_member m
                      JOIN compound_set s ON s.set_id = m.set_id
                     WHERE {' AND '.join(member)})"""
            )
        if family is not None:
            # measured against, not active against. deciding what counts as
            # active needs a threshold and a unit, and the database picks
            # neither, so a counter screen at >30 uM is still an answer to
            # "what is available for this family"
            where.append(
                f"""c.inchikey IN (
                    SELECT g.inchikey FROM bioactivity_group g
                      JOIN target_uniprot tu ON tu.target_id = g.target_id
                      JOIN uniprot u ON u.uniprot_id = tu.uniprot_id
                     WHERE {FAMILY_MATCH})"""
            )
            params.append(family)

        return self.read(
            f"""
            SELECT c.inchikey, c.name, c.smiles,
                   COUNT(s.set_id) AS n_sets,
                   GROUP_CONCAT(s.name, ', ') AS sets,
                   (SELECT COUNT(DISTINCT g.target_id) FROM bioactivity_group g
                     WHERE g.inchikey = c.inchikey) AS n_targets
              FROM compound c
              LEFT JOIN compound_set_member cm ON cm.inchikey = c.inchikey
              LEFT JOIN compound_set s ON s.set_id = cm.set_id
              {"WHERE " + " AND ".join(where) if where else ""}
             GROUP BY c.inchikey
             ORDER BY n_targets DESC, c.name
        """,
            *params,
        )

    def families(self, like=None):
        """Every level of the UniProt classification, with what sits under it."""
        # the stored value is the whole hierarchy, outermost first, so
        # "Small GTPase superfamily, Ras family" is two answers and not one.
        # splitting it here means you can ask for either level by name
        rows = self.read(
            """
            SELECT u.superfamily,
                   COUNT(DISTINCT u.uniprot_id) AS proteins,
                   COUNT(DISTINCT tu.target_id) AS targets
              FROM uniprot u
              LEFT JOIN target_uniprot tu ON tu.uniprot_id = u.uniprot_id
             WHERE u.superfamily IS NOT NULL AND u.superfamily != ''
             GROUP BY u.superfamily
        """
        )
        if rows.empty:
            return rows.assign(family=None)
        out = (
            rows.assign(family=rows["superfamily"].str.split(", "))
            .explode("family")
            .groupby("family", as_index=False)
            .agg(proteins=("proteins", "sum"), targets=("targets", "sum"))
            .sort_values(["proteins", "family"], ascending=[False, True])
        )
        if like is not None:
            # on whole words, or asking for "Ras" also finds "transferase"
            out = out[
                out["family"].str.contains(
                    rf"\b{re.escape(like)}\b", case=False, regex=True
                )
            ]
        return out.reset_index(drop=True)

    def sets(self):
        return self.read(
            """
            SELECT s.set_id, s.name, s.category, s.source_db,
                   COUNT(m.inchikey) AS compounds
              FROM compound_set s
              LEFT JOIN compound_set_member m ON m.set_id = s.set_id
             GROUP BY s.set_id
             ORDER BY compounds DESC, s.name
        """
        )

    def compound_sets(self, key):
        # the other direction: given a compound, which collections claim it
        return self.read(
            "SELECT s.set_id, s.name, s.category, s.source_db, s.description "
            "  FROM compound_set_member m "
            "  JOIN compound_set s ON s.set_id = m.set_id "
            " WHERE m.inchikey = ? ORDER BY s.category, s.name",
            self.compound_key(key),
        )

    def targets(self, key=None, family=None):
        # target and uniprot share no column, target_uniprot is the link.
        # one row per target and accession, so a complex comes back as
        # several rows with the same target_id
        where, params = [], []
        if key is not None:
            ids = self.targets_for(key)
            where.append(f"target_id IN ({','.join('?' * len(ids))})")
            params += ids
        if family is not None:
            # a complex is in the family if any of its members is, and it still
            # comes back with all its members, not only the classified one
            where.append(
                f"""target_id IN (
                    SELECT tu.target_id FROM target_uniprot tu
                      JOIN uniprot u ON u.uniprot_id = tu.uniprot_id
                     WHERE {FAMILY_MATCH})"""
            )
            params.append(family)
        sql = "SELECT * FROM target_flat"
        if where:
            sql += " WHERE " + " AND ".join(where)
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

    def bioactivities(self, compound=None, target=None, family=None):
        # the view is the join; this only filters it
        where, params = [], []
        if compound is not None:
            where.append("inchikey = ?")
            params.append(self.compound_key(compound))
        if target is not None:
            ids = self.targets_for(target)
            where.append(f"target_id IN ({','.join('?' * len(ids))})")
            params += ids
        if family is not None:
            where.append(
                f"""target_id IN (
                    SELECT tu.target_id FROM target_uniprot tu
                      JOIN uniprot u ON u.uniprot_id = tu.uniprot_id
                     WHERE {FAMILY_MATCH})"""
            )
            params.append(family)
        sql = "SELECT * FROM bioactivity_flat"
        if where:
            sql += " WHERE " + " AND ".join(where)
        return self.read(sql + " ORDER BY id", *params)

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
        key = str(key).strip()
        rows = self.conn.execute(
            "SELECT DISTINCT tu.target_id FROM target_uniprot tu "
            "  JOIN uniprot u ON u.uniprot_id = tu.uniprot_id "
            " WHERE tu.uniprot_id = ? OR u.hgnc = ? ORDER BY tu.target_id",
            (key.upper(), key.upper()),
        ).fetchall()
        # and by name as well, not instead. not every uniprot row carries an
        # HGNC symbol, and 302 protein targets have no accession at all, so
        # matching only on the accession side silently drops them: asking for
        # EGFR by symbol finds a 152 row target and misses the 3263 row one
        named = self.conn.execute(
            "SELECT target_id FROM target WHERE name = ? COLLATE NOCASE",
            (key,),
        ).fetchall()
        return sorted({r[0] for r in rows} | {r[0] for r in named})

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

    def add_set(self, name, category, source_db=None, description=None):
        found = self.one("SELECT set_id FROM compound_set WHERE name = ?", name)
        if found is not None:
            return found
        return self.insert(
            "compound_set",
            name=name,
            category=category,
            source_db=source_db or None,
            description=description or None,
        )

    def add_set_member(self, set_id, inchikey):
        self.conn.execute(
            "INSERT OR IGNORE INTO compound_set_member (set_id, inchikey) "
            "VALUES (?, ?)",
            (set_id, inchikey.strip().upper()),
        )
        return set_id

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

    def set_superfamily(self, uniprot_id, superfamily):
        uniprot_id = uniprot_id.strip().upper()
        self.conn.execute(
            "UPDATE uniprot SET superfamily = ? WHERE uniprot_id = ?",
            (superfamily or None, uniprot_id),
        )

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
