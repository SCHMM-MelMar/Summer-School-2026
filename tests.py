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

# every datapoint names the database it came from, and where the source gives
# a resolvable id the api hands back a link instead of leaving it to be guessed
assert rows["source_db"].notna().all()
assert set(rows["source_db"]) == {"literature", "opnMe", "Probes & Drugs", "in-house"}
lit = rows[rows.source_db == "literature"]
assert lit["source"].str.startswith("PMID:").all()
assert lit["source_url"].str.startswith("https://doi.org/10.").all()
assert rows[rows.source_db == "in-house"]["source_url"].isna().all()
assert db.sources()["measurements"].sum() == len(rows)

# one directory is one set, named after itself, so every compound can say where
# it came from. this is one row per compound, not one row per membership
compounds = db.compounds()
assert len(compounds) == counts["compound"]
assert len(db.table("compound_flat")) == counts["compound_set_member"]
assert db.sets()["name"].tolist() == ["template"]
assert db.sets()["compounds"].tolist() == [3]
assert len(db.compounds(set="template")) == 3
assert len(db.compounds(set="not a set")) == 0
assert compounds.set_index("inchikey").loc[JQ1, "sets"] == "template"

# the other direction, from a compound to what claims it
assert db.compound_sets("BI-2536")["name"].tolist() == ["template"]

# the UniProt classification is a property of one protein, and it is stored as
# the whole hierarchy, so either level answers by name
db.set_superfamily("P53350", "Protein kinase superfamily, Ser/Thr protein kinase family")
assert db.targets(family="Ser/Thr protein kinase family").target_id.tolist() == [1]
assert db.targets(family="Protein kinase superfamily").target_id.tolist() == [1]
assert len(db.targets(family="Ras family")) == 0
assert db.families()["family"].tolist() == [
    "Protein kinase superfamily",
    "Ser/Thr protein kinase family",
]
# on whole words, so asking for Ras does not find transferase
assert len(db.families(like="kinase")) == 2
assert len(db.families(like="Ras")) == 0
# measured against, not active against: Olaparib is in because the template
# counter-screens it on PLK1 at Kd > 30 uM, which is a real answer to
# "what is available", just not a hit
assert sorted(db.compounds(family="Protein kinase superfamily").name) == [
    "BI-2536",
    "Olaparib",
]

# target and uniprot are two tables with no shared column, targets() joins them
flat = db.targets()
assert set(flat["target_id"]) == set(db.table("target")["target_id"])
assert set(db.targets("PARP1")["uniprot_id"]) == {"P09874", "Q9UGN5", "Q9Y6F1"}
assert (db.targets("P06493")["type"] == "complex").all()

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

# structural similarity. RDKit is optional, so this is skipped rather than
# failed when it is not installed
from probedb import similarity

try:
    similarity.fingerprints(db)
except ImportError as missing:
    logger.warning(f"skipping the similarity checks: {missing}")
else:
    assert len(similarity.fingerprints(db)) == 3
    scored = similarity.pairs(db)
    assert len(scored) == 3  # one row per pair, never the mirror image
    assert (scored.tanimoto < 0.5).all()  # three unrelated chemotypes
    assert not scored.same_skeleton.any()
    assert similarity.neighbours(db, "BI-2536", n=1).tanimoto.iloc[0] < 0.5
    square = similarity.matrix(db)
    assert square.shape == (3, 3)
    assert (square.to_numpy().diagonal() == 1.0).all()

    # the two decisions this module makes, checked rather than asserted in a
    # comment. counts, because a binary fingerprint cannot count carbons:
    bench = ProbeDB(":memory:", create=True)
    bench.add_compound("POULHZVOKOAJMA-UHFFFAOYSA-N", "CCCCCCCCCCCC(=O)O", "lauric acid")
    bench.add_compound(
        "UKMSUNONTOPOIO-UHFFFAOYSA-N", "CCCCCCCCCCCCCCCCCCCCCC(=O)O", "behenic acid"
    )
    homologues = similarity.pairs(bench).tanimoto.iloc[0]
    assert homologues < 0.7, f"C12 and C22 scored {homologues}, counts are off"

    # and chirality, because a probe and its inactive enantiomer are the pair
    # that most needs telling apart
    bench = ProbeDB(":memory:", create=True)
    JQ1_MINUS = "DNVXATUJJDPFDM-QGZVFWFLSA-N"
    bench.add_compound(JQ1, db.one("SELECT smiles FROM compound WHERE inchikey = ?", JQ1))
    bench.add_compound(
        JQ1_MINUS,
        "Cc1sc2c(c1C)C(c1ccc(Cl)cc1)=N[C@H](CC(=O)OC(C)(C)C)c1nnc(C)n1-2",
    )
    enantiomers = similarity.pairs(bench).iloc[0]
    assert enantiomers.tanimoto < 1.0, "the JQ1 enantiomers scored identical"
    assert enantiomers.same_skeleton  # same constitution, different stereochemistry

sources = [
    d for d in sorted(STAGING.iterdir()) if d.is_dir() and not d.name.startswith("_")
]
if not sources:
    logger.warning("no source directories in staging/, checked the template only")
else:
    # staging/ is where everybody drops their source, so a directory that is
    # still being worked on is reported and skipped rather than failing the run
    full, loaded = ProbeDB(":memory:", create=True), []
    for directory in sources:
        try:
            load(full, directory, source=directory.name)
            loaded.append(directory.name)
        except ValueError as problem:
            logger.warning(f"skipped {directory.name}\n{problem}")
    logger.info(f"loaded {loaded}: {dict(full.counts().values)}")

logger.info("ok")
