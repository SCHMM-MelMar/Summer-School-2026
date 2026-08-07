from pathlib import Path

from loguru import logger

from probedb import ProbeDB
from loader import load_all, load_families

STAGING = Path(__file__).resolve().parent.parent / "staging"

db = ProbeDB("probe.db", create=True)

reports = load_all(db, STAGING)

if not reports:
    logger.warning(
        f"no source directories in {STAGING} -- run `roche-css staging`, "
        f"or copy staging/_template and fill it in"
    )
    raise SystemExit(1)

for report in reports:
    counts = {
        k: report[k] for k in ("compounds", "targets", "complexes", "bioactivities")
    }
    logger.info(f"{report['directory']}  {counts}")

updated = load_families(db)
logger.info(f"uniprot.superfamily backfilled for {updated} accessions")

logger.info("\n" + db.counts().to_string(index=False))
