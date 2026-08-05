from pathlib import Path

from loguru import logger

from probedb import ProbeDB
from loader import load, load_all, validate

STAGING = Path(__file__).resolve().parent / "staging"

BI_2536 = "XQVVPGYIWAGRNI-JOCHJYFZSA-N"
JQ1 = "DNVXATUJJDPFDM-KRWDZBQOSA-N"
PARP = ["P09874", "Q9UGN5", "Q9Y6F1"]

assert validate(STAGING / "_template") == [], validate(STAGING / "_template")

# a ChEMBL id belongs to one compound, and validate says so before a load starts
import csv, shutil, tempfile
scratch = Path(tempfile.mkdtemp())
shutil.copytree(STAGING / "_template", scratch / "d")
rows = list(csv.DictReader(open(scratch / "d" / "compound.tsv"), delimiter="\t"))
rows[1]["chembl_id"] = rows[0]["chembl_id"]
with open(scratch / "d" / "compound.tsv", "w", newline="") as fh:
    writer = csv.DictWriter(fh, rows[0].keys(), delimiter="\t", lineterminator="\n")
    writer.writeheader()
    writer.writerows(rows)
assert any("belongs to one compound" in p for p in validate(scratch / "d"))
shutil.rmtree(scratch)

db = ProbeDB(":memory:", create=True)
load(db, STAGING / "_template", source="template")
rows = db.bioactivities()

counts = dict(db.counts().values)
assert counts["compound"] == 3
assert counts["target"] == 6
assert counts["bioactivity"] == 13

# a target is one row whatever it is made of
assert sorted(db.read("SELECT DISTINCT type FROM target")["type"]) == [
    "complex",
    "family",
    "protein",
]

# a member and the group it belongs to are different targets
assert len(db.targets_for("P09874")) == 2
family = db.find_target("family", PARP)
assert db.find_target("family", reversed(PARP)) == family
assert db.find_target("family", PARP[:2]) is None
assert db.one("SELECT COUNT(*) FROM target_flat WHERE target_id = ?", family) == 3

# compounds are searchable by name
assert db.find("JQ1").inchikey.iloc[0] == JQ1
assert db.compound_key("(+)-JQ1") == JQ1
assert db.compound_key("CHEMBL1957266") == JQ1  # by ChEMBL id

# chembl maps an InChIKey to the ChEMBL id itself
assert db.one("SELECT chembl_id FROM chembl WHERE inchikey = ?", JQ1) == "CHEMBL1957266"
assert db.one("SELECT COUNT(*) FROM chembl") == 3

# how a number was measured sits on the measurement itself
assert set(rows.assay_type.dropna()) == {"biochemical", "cell", "binding"}
assert rows[rows.cell_line == "HeLa"].bioactivity_type.iloc[0] == "EC50"

# percent readouts are not a bioactivity type and are not loaded anywhere yet
assert "%" not in set(rows["unit"])
assert not any(str(t).startswith("%") for t in rows["bioactivity_type"])

# units are kept as reported, nothing is converted anywhere
assert set(rows["unit"]) == {"nM", "uM", "-log(M)"}
assert {"IC50", "pIC50"} <= set(rows["bioactivity_type"])

# censored values are stored rather than dropped
assert (rows["relation"] != "=").sum() == 2

# an unknown mode of action is empty, not null, so the group foreign key holds
assert db.one("SELECT COUNT(*) FROM bioactivity WHERE moa IS NULL") == 0
assert db.one("SELECT COUNT(*) FROM bioactivity WHERE moa = ''") == 2
assert (
    db.one(
        "SELECT COUNT(*) FROM bioactivity b WHERE NOT EXISTS "
        "(SELECT 1 FROM bioactivity_group g WHERE g.inchikey = b.inchikey "
        " AND g.target_id = b.target_id AND g.moa = b.moa)"
    )
    == 0
)

# repeated readings stay separate rows
brd2 = rows[rows.target == "Bromodomain-containing protein 2"]
assert (brd2["source"] == "PMID:35007061").sum() == 2
assert set(brd2["bioactivity_type"]) == {"IC50", "Kd"}

# one pair, three source databases
plk1 = rows[
    (rows.inchikey == BI_2536) & (rows.target_id == db.targets_for("P53350")[0])
]
assert len(set(plk1["source_db"])) == 3 and len(plk1) == 4

# every datapoint names the database it came from. A record id is only there
# when the source has one: a number read from a paper is identified by the PMID.
assert rows["source_db"].notna().all()
assert set(rows["source_db"]) == {"literature", "opnMe", "Probes & Drugs", "in-house"}
assert rows[rows.source_db == "literature"]["source"].str.startswith("PMID:").all()
assert rows[rows.source_db == "literature"]["source_xref"].isna().all()

# a target type outside the vocabulary is refused
try:
    db.add_target("nonsense", "x", [("P00000", None, None)])
    raise AssertionError("bad target type was accepted")
except Exception:
    pass

try:
    db.add_compound("not-an-inchikey")
    raise AssertionError("malformed InChIKey was accepted")
except ValueError:
    pass

# loading the same directory twice adds nothing: an exact repeat of a
# measurement is a duplicate, not a second piece of evidence
report = load(db, STAGING / "_template", source="template")
again = dict(db.counts().values)
assert again == counts
assert report["bioactivities"] == 0
assert report["duplicates_skipped"] == counts["bioactivity"]

sources = [
    d for d in sorted(STAGING.iterdir()) if d.is_dir() and not d.name.startswith("_")
]
if not sources:
    logger.warning("no source directories in staging/, checked the template only")
else:
    full = ProbeDB(":memory:", create=True)
    load_all(full, STAGING)
    logger.info(f"loaded {len(sources)} source(s): " f"{dict(full.counts().values)}")

logger.info("ok")
