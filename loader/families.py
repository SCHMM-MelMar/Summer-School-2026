import csv
from pathlib import Path

REFERENCE = (
    Path(__file__).resolve().parent.parent / "reference" / "uniprot_protein_families.tsv"
)


def load_families(db, path=REFERENCE):
    """Backfill uniprot.superfamily from the UniProt protein-family reference
    file. Only touches accessions already present in the uniprot table."""
    path = Path(path)
    if not path.exists():
        return 0

    known = {r[0] for r in db.conn.execute("SELECT uniprot_id FROM uniprot")}

    updated = 0
    with open(path, newline="") as handle:
        for row in csv.DictReader(handle, delimiter="\t"):
            uniprot_id = row["uniprot_id"].strip().upper()
            superfamily = row.get("protein_family", "").strip()
            if uniprot_id not in known or not superfamily:
                continue
            db.set_superfamily(uniprot_id, superfamily)
            updated += 1

    db.commit()
    return updated
