"""Extract compound/target/uniprot/bioactivity TSVs from the Probes & Drugs
sqlite dump, scoped to compoundset 395 (High-quality chemical probes / HQCP).
See gap_analysis.md for the design decisions this implements.

Disclaimer: generated with Claude Code (v2.1.179, model claude-sonnet-5),
with parts reviewed and cleaned up manually.

Usage:
    python3 extract.py [--db pd_02_2025_dump.sqlite] [--out output]
"""

import argparse
import csv
import sqlite3
from collections import defaultdict
from itertools import zip_longest
from pathlib import Path

# P&D's own curated compound sets, defined by the compoundset table and
# populated via compoundtocompoundset. 395 = "High-quality chemical probes"
COMPOUNDSET_IDS = (395,)  # high-quality chemical probes (HQCP)
CHEMBL_EXTERNALDB_ID = 2

TYPE_MAP = {
    1: "family",     # protein complex group
    2: "complex",    # protein complex
    3: "family",     # selectivity group (counter-screening panels, not obligate complexes)
    4: "family",     # protein family
    5: "protein",    # single protein
    6: "ppi",        # protein-protein interaction
    12: "organism",  # whole-organism phenotypic assay -- no uniprot accession exists
    13: "cell_line",  # whole-cell-line phenotypic assay -- no uniprot accession exists
    14: "protein",   
    }


CELLLINE_COLLISIONS = {"B-95-8 cell line", "BSC-1", "CHO", "NG108-15", "TPC1"}


RELATION_MAP = {
    "=": "=", ">": ">", "<": "<", ">=": ">=", "<=": "<=",
    "median": "~", "max": "~", "min": "~",
    "-": "", "": "",
}


def connect(db_path):
    """Open the dump and materialize the compound-id scope once as a temp
    table, so every later query just joins against `scope` instead of
    repeating the compoundtocompoundset filter (and risking it drifting out
    of sync between queries)."""
    conn = sqlite3.connect(db_path)
    conn.row_factory = sqlite3.Row
    conn.execute(
        "CREATE TEMP TABLE scope AS "
        "SELECT DISTINCT compound_id FROM compoundtocompoundset "
        f"WHERE compoundset_id IN ({','.join('?' * len(COMPOUNDSET_IDS))})",
        COMPOUNDSET_IDS,
    )
    conn.execute("CREATE INDEX idx_scope_compound_id ON scope (compound_id)")
    return conn


def extract_compounds(conn):
    """Returns (compound_id -> inchikey) for valid compounds, and writes rows
    for compound.tsv (one per distinct inchikey)."""

    chembl_ids = defaultdict(dict)
    for row in conn.execute(
        "SELECT compound_id, ligand_id FROM compoundtoexternaldb "
        f"WHERE externaldb_id = {CHEMBL_EXTERNALDB_ID} "
        "AND compound_id IN (SELECT compound_id FROM scope)"
    ):
        chembl_ids[row["compound_id"]][row["ligand_id"]] = None

    compound_id_to_inchikey = {}
    seen_inchikeys = {}
    rows = []
    skipped_dupe = 0

    for row in conn.execute(
        "SELECT c.compoundid, c.inchikey, c.smiles, c.name "
        "FROM compound c JOIN scope s ON s.compound_id = c.compoundid "
        "ORDER BY c.compoundid"
    ):
        # no format check -- every compound in this scope is expected to
        # already carry a well-formed inchikey by design
        inchikey = (row["inchikey"] or "").upper()

        compound_id_to_inchikey[row["compoundid"]] = inchikey
        if inchikey in seen_inchikeys:
            skipped_dupe += 1
            continue
        seen_inchikeys[inchikey] = True
        rows.append(
            {
                "inchikey": inchikey,
                "smiles": row["smiles"] or "",
                "chembl_id": ";".join(chembl_ids.get(row["compoundid"], {})),
                "name": row["name"] or "",
            }
        )

    print(f"compound.tsv: {len(rows)} rows ({skipped_dupe} duplicate inchikey collapsed)")
    return compound_id_to_inchikey, rows


def extract_targets_and_uniprot(conn):
    """Returns (target_id -> target_key) for targets actually referenced by
    scoped activity rows, and writes rows for target.tsv and uniprot.tsv."""
    target_rows, uniprot_rows = [], []
    declared = {}
    target_id_to_key = {}
    skipped_unmapped_type = 0

    # Only pull targets that scoped activity actually points at (via the
    # JOIN on activity), not every target row in the DB. P&D's `target`
    # table has ~37k rows total, the vast majority irrelevant to our 922
    # compounds.
    query = (
        "SELECT DISTINCT t.targetid, t.name, t.targettype_id, t.uniprotid, "
        "       t.gene_name, o.name AS species "
        "FROM target t "
        "JOIN activity a ON a.target_id = t.targetid "
        "JOIN scope s ON s.compound_id = a.compound_id "
        "LEFT JOIN organism o ON o.organismid = t.organism_id"
    )
    for row in conn.execute(query):
        target_type = TYPE_MAP.get(row["targettype_id"])
        if target_type is None:
            skipped_unmapped_type += 1
            continue

        name = row["name"] or ""
        target_key = name
        if target_type == "cell_line" and name in CELLLINE_COLLISIONS:
            target_key = f"cellline:{name}"
        target_id_to_key[row["targetid"]] = target_key

        prior = declared.get(target_key)
        if prior is None:
            declared[target_key] = (target_type, name)
            target_rows.append({"target_key": target_key, "type": target_type, "name": name})
        elif prior != (target_type, name):
            print(f"  ! target_key {target_key!r} redeclared with a different type/name, kept first")

        if target_type in ("cell_line", "organism"):
            continue  # no accession -- phenotypic target, zero uniprot.tsv rows

        accessions = [a for a in (row["uniprotid"] or "").split(",") if a]
        genes = [g for g in (row["gene_name"] or "").split(",") if g]
        if len(accessions) != len(genes):
            print(f"  ! target {name!r}: {len(accessions)} accessions vs {len(genes)} gene names, hgnc blank past the shorter list")

        for accession, gene in zip_longest(accessions, genes, fillvalue=""):
            if not accession:
                continue  
            uniprot_rows.append(
                {
                    "uniprot_id": accession,
                    "target_key": target_key,
                    "hgnc": gene,
                    "species": row["species"] or "",
                }
            )

    # Used later in extract_bioactivity to recognize a cell_line-type target
    # by its target_key alone, without re-querying the DB.
    target_key_to_type = {key: type_ for key, (type_, _name) in declared.items()}

    print(
        f"target.tsv: {len(target_rows)} rows ({skipped_unmapped_type} activity-linked "
        f"targets skipped, unmapped targettype)"
    )
    print(f"uniprot.tsv: {len(uniprot_rows)} rows")
    return target_id_to_key, target_key_to_type, target_rows, uniprot_rows


def resolve_moa(conn):
    """(compound_id, target_id) -> actiontype_id, per the join rule in
    gap_analysis.md item 4: prefer primary_target=1, else agreement, else
    lowest compoundactionid; log genuine conflicts.

    MOA (mechanism of action, e.g. "inhibitor") is recorded in P&D as a fact
    about a (compound, target) *pair* -- compoundaction -- completely
    separate from the activity table's individual measurements. There's no
    moa column on `activity` itself, so this has to be a standalone lookup
    joined in afterwards. A pair can have more than one compoundaction row
    (recorded by different curated source collections), occasionally
    disagreeing on the actiontype -- e.g. one source calls it "inhibitor",
    another "antagonist" for the same compound+target. Those conflicts are
    real and rare (checked: <1% of pairs in this scope), so we pick
    deterministically rather than silently overwriting with whichever row the
    SQL engine happens to return last.
    """
    pairs = defaultdict(list)
    for row in conn.execute(
        "SELECT compoundactionid, compound_id, target_id, actiontype_id, primary_target "
        "FROM compoundaction WHERE compound_id IN (SELECT compound_id FROM scope)"
    ):
        pairs[(row["compound_id"], row["target_id"])].append(row)

    moa = {}
    conflicts = 0
    for key, rows in pairs.items():
        primary = [r for r in rows if r["primary_target"] == 1 and r["actiontype_id"]]
        if primary:
            moa[key] = primary[0]["actiontype_id"]
            continue
        types = {r["actiontype_id"] for r in rows if r["actiontype_id"]}
        if len(types) <= 1:
            moa[key] = next(iter(types)) if types else ""
        else:
            conflicts += 1
            moa[key] = sorted(rows, key=lambda r: r["compoundactionid"])[0]["actiontype_id"]

    print(f"MOA lookup: {len(moa)} (compound, target) pairs, {conflicts} conflicting resolved to lowest id")
    return moa


def reference_label(ref):
    """Turn one `reference` row into the human-readable string that goes in
    bioactivity.tsv's `source` column. Preference order follows what's
    actually resolvable/citable, falling back progressively for the
    referencetypes that don't carry a PMID/DOI/patent (dataset, webpage,
    book, pubchem_bioassay)."""
    if ref["pubmed_id"]:
        return f"PMID:{ref['pubmed_id']}"
    if ref["doi"]:
        return f"DOI:{ref['doi']}"
    if ref["patent_id"]:
        return f"Patent:{ref['patent_id']}"
    if ref["url"]:
        return ref["url"]
    if ref["title"]:
        return ref["title"][:120]
    return f"reference:{ref['referenceid']}"


def extract_bioactivity(conn, compound_id_to_inchikey, target_id_to_key, target_key_to_type, moa):
    # Preload every (activity_id -> [reference rows]) once, rather than
    # querying per-activity inside the main loop below
    references = defaultdict(list)
    ref_query = (
        "SELECT atr.activity_id, r.referenceid, r.pubmed_id, r.doi, r.patent_id, r.url, r.title "
        "FROM activitytoreference atr "
        "JOIN reference r ON r.referenceid = atr.reference_id "
        "JOIN activity a ON a.activityid = atr.activity_id "
        "JOIN scope s ON s.compound_id = a.compound_id"
    )
    for row in conn.execute(ref_query):
        references[row["activity_id"]].append(row)

    cell_lines = {}
    for row in conn.execute("SELECT celllineid, name FROM cell_line"):
        cell_lines[row["celllineid"]] = row["name"]

    rows = []
    skipped_compound, skipped_target = 0, 0
    unmapped_relations = defaultdict(int)

    query = (
        "SELECT a.activityid, a.compound_id, a.target_id, a.activity_type, "
        "       a.activity_value, a.value_type, a.cellline_id, a.cellline_custom, "
        "       a.activity_unit_id "
        "FROM activity a JOIN scope s ON s.compound_id = a.compound_id"
    )
    for row in conn.execute(query):
        inchikey = compound_id_to_inchikey.get(row["compound_id"])
        if inchikey is None:
            # defensive only -- compound_id_to_inchikey covers every scoped
            # compound unconditionally, so this shouldn't be reachable
            skipped_compound += 1
            continue
        target_key = target_id_to_key.get(row["target_id"])
        if target_key is None:
            skipped_target += 1
            continue

        relation = row["value_type"] if row["value_type"] is not None else ""
        if relation not in RELATION_MAP:
            unmapped_relations[relation] += 1
            relation = ""
        else:
            relation = RELATION_MAP[relation]
        cell_line = cell_lines.get(row["cellline_id"], "") or row["cellline_custom"] or ""
        if not cell_line and target_key_to_type.get(target_key) == "cell_line":
            cell_line = target_key.removeprefix("cellline:")

        base = {
            "inchikey": inchikey,
            "target_key": target_key,
            "moa": moa.get((row["compound_id"], row["target_id"]), ""),
            "bioactivity_type": row["activity_type"] or "",
            "relation": relation,
            "value": "" if row["activity_value"] is None else row["activity_value"],
            "unit": row["activity_unit_id"] or "",
            "assay_type": "",
            "assay_description": "",
            "cell_line": cell_line,
            "concentration": "",
            "concentration_unit": "",
            "source_db": "Probes & Drugs",
        }

        refs = references.get(row["activityid"], [])
        if not refs:
            rows.append({**base, "source": "", "source_xref": "", "xref_id": ""})
        else:
            for ref in refs:
                rows.append(
                    {
                        **base,
                        "source": reference_label(ref),
                        # source_xref is P&D's own internal activity id
                        "source_xref": f"activity:{row['activityid']}",
                        "xref_id": "",
                    }
                )

    if unmapped_relations:
        print(f"  ! unmapped value_type values (dropped to blank relation): {dict(unmapped_relations)}")
    print(
        f"bioactivity.tsv: {len(rows)} rows "
        f"({skipped_compound} skipped, compound not in scope; "
        f"{skipped_target} skipped, target type unmapped)"
    )
    return rows


def write_tsv(path, columns, rows):
    path.parent.mkdir(parents=True, exist_ok=True)
    with open(path, "w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=columns, delimiter="\t")
        writer.writeheader()
        writer.writerows(rows)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--db", default=Path(__file__).parent / "pd_02_2025_dump.sqlite")
    parser.add_argument("--out", default=Path(__file__).parent / "output")
    args = parser.parse_args()

    out = Path(args.out)
    conn = connect(args.db)

    compound_id_to_inchikey, compound_rows = extract_compounds(conn)
    target_id_to_key, target_key_to_type, target_rows, uniprot_rows = extract_targets_and_uniprot(conn)
    moa = resolve_moa(conn)
    bioactivity_rows = extract_bioactivity(
        conn, compound_id_to_inchikey, target_id_to_key, target_key_to_type, moa
    )

    write_tsv(out / "compound.tsv", ["inchikey", "smiles", "chembl_id", "name"], compound_rows)
    write_tsv(out / "target.tsv", ["target_key", "type", "name"], target_rows)
    write_tsv(out / "uniprot.tsv", ["uniprot_id", "target_key", "hgnc", "species"], uniprot_rows)
    write_tsv(
        out / "bioactivity.tsv",
        [
            "inchikey", "target_key", "moa", "bioactivity_type", "relation", "value", "unit",
            "assay_type", "assay_description", "cell_line", "concentration", "concentration_unit",
            "source_db", "source", "source_xref", "xref_id",
        ],
        bioactivity_rows,
    )
    conn.close()
    print(f"done -> {out}")


if __name__ == "__main__":
    main()
