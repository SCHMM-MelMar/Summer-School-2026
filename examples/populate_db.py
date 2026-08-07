import argparse
from pathlib import Path

from loguru import logger

from probedb import ProbeDB
from loader import load, load_families

STAGING = Path(__file__).resolve().parent.parent / "staging"

LOADED = ("compounds", "targets", "complexes", "sets", "bioactivities")
SKIPPED = (
    "targets_skipped",
    "compounds_skipped",
    "accessions_skipped",
    "bioactivities_skipped",
    "duplicates_skipped",
)

parser = argparse.ArgumentParser(description="build probe.db from staging/")
parser.add_argument("--out", default="probe.db", help="where to write the database")
parser.add_argument(
    "--lenient",
    action="store_true",
    help="load what each directory can represent instead of refusing the whole "
    "directory over rows the schema has no home for. everything dropped is "
    "still counted and reported",
)
args = parser.parse_args()

path = Path(args.out)
if path.exists():
    logger.warning(f"{path} already exists, delete it first for a clean rebuild")
    raise SystemExit(1)

db = ProbeDB(path, create=True)

sources = [
    d for d in sorted(STAGING.iterdir()) if d.is_dir() and not d.name.startswith("_")
]
if not sources:
    logger.warning(
        f"no source directories in {STAGING} -- run `roche-css staging`, "
        f"or copy staging/_template and fill it in"
    )
    raise SystemExit(1)

# staging/ is where everybody drops their source, so a directory that does not
# validate is reported and skipped instead of taking the whole run down
for directory in sources:
    try:
        report = load(db, directory, source=directory.name, strict=not args.lenient)
    except ValueError as problem:
        logger.warning(f"skipped {directory.name}\n{problem}")
        continue
    loaded = {k: report[k] for k in LOADED}
    logger.info(f"{directory.name}  {loaded}")
    dropped = {k: report[k] for k in SKIPPED if report[k]}
    if dropped:
        logger.warning(f"{directory.name}  dropped {dropped}")

# the family classification is the same for everybody, so it comes from one
# reference file rather than from whichever source happened to mention it
updated = load_families(db)
logger.info(f"uniprot.superfamily filled in for {updated} accessions")

logger.info("\n" + db.counts().to_string(index=False))
logger.info(f"wrote {path.resolve()}")
