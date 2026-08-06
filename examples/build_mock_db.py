from pathlib import Path

from loguru import logger

from probedb import ProbeDB
from loader import load
from loader.load import STAGING

DB_PATH = Path(__file__).resolve().parent.parent / "probe.db"

DB_PATH.unlink(missing_ok=True)

db = ProbeDB(DB_PATH, create=True)

report = load(db, STAGING / "_template", source="template")
db.commit()

logger.info(report)
logger.info("\n" + db.counts().to_string(index=False))
logger.info(f"wrote {DB_PATH}")
