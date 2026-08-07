"""Checks on a built database.

Loading reports what each directory contributed. This asks the opposite
question: now that six sources are in one place, did they actually merge, or
are they six databases sharing a file? Run it after populate_db.py.

    python examples/qc_db.py probe.db
"""

import argparse
import re

from loguru import logger

from probedb import ProbeDB
from probedb.db import INCHIKEY
from probedb.schema import vocabulary

UNIPROT = re.compile(
    r"^([OPQ][0-9][A-Z0-9]{3}[0-9]|[A-NR-Z][0-9]([A-Z][A-Z0-9]{2}[0-9]){1,2})$"
)

parser = argparse.ArgumentParser(description=__doc__)
parser.add_argument("path", nargs="?", default="probe.db")
args = parser.parse_args()

CONCENTRATION_UNITS = ("M", "mM", "uM", "µM", "nM", "pM", "fM")

db = ProbeDB(args.path)
checks = []
notes = []


def check(ok, what, detail=""):
    """Something that would mean the database is wrong."""
    checks.append((bool(ok), what, detail))


def note(what, detail=""):
    """Something worth knowing that is the source data's business, not ours."""
    notes.append((what, detail))


def count(sql, *params):
    return db.one(sql, *params)


# -- nothing is dangling ---------------------------------------------------

broken = list(db.conn.execute("PRAGMA foreign_key_check"))
check(not broken, "every foreign key resolves", f"{len(broken)} broken")

check(
    count("SELECT COUNT(*) FROM compound WHERE inchikey IS NULL OR inchikey = ''") == 0,
    "every compound has a key",
)
bad_keys = [
    k
    for (k,) in db.conn.execute("SELECT inchikey FROM compound")
    if not INCHIKEY.match(k)
]
check(not bad_keys, "every compound key is a well-formed InChIKey", str(bad_keys[:5]))

bad_acc = [
    a
    for (a,) in db.conn.execute("SELECT uniprot_id FROM uniprot")
    if not UNIPROT.match(a)
]
check(
    not bad_acc,
    "every accession is a well-formed UniProt id",
    f"{len(bad_acc)}: {bad_acc[:5]}",
)

# a target with neither a name nor an accession cannot be identified again
nameless = count(
    "SELECT COUNT(*) FROM target t WHERE (t.name IS NULL OR t.name = '') "
    "AND NOT EXISTS (SELECT 1 FROM target_uniprot tu WHERE tu.target_id = t.target_id)"
)
check(
    nameless == 0, "every target has a name or an accession", f"{nameless} have neither"
)

# the group table is what an aggregation would run over, so a measurement
# outside it would silently disappear from every summary
orphan_groups = count(
    "SELECT COUNT(*) FROM bioactivity b WHERE NOT EXISTS ("
    "  SELECT 1 FROM bioactivity_group g WHERE g.inchikey = b.inchikey"
    "   AND g.target_id = b.target_id AND g.moa = b.moa)"
)
check(orphan_groups == 0, "every measurement is in a group", f"{orphan_groups} are not")

check(
    count("SELECT COUNT(*) FROM bioactivity WHERE moa IS NULL") == 0,
    "an unknown mode of action is empty and not null",
)

# -- vocabularies held ------------------------------------------------------

for column, table in (("type", "target"), ("relation", "bioactivity")):
    allowed = vocabulary(column)
    marks = ",".join("?" * len(allowed))
    off = count(
        f"SELECT COUNT(*) FROM {table} WHERE {column} IS NOT NULL "
        f"AND {column} NOT IN ({marks})",
        *sorted(allowed),
    )
    check(off == 0, f"{table}.{column} stays inside its vocabulary", f"{off} outside")

off = count(
    "SELECT COUNT(*) FROM compound_set WHERE category NOT IN ({})".format(
        ",".join("?" * len(vocabulary("category")))
    ),
    *sorted(vocabulary("category")),
)
check(off == 0, "compound_set.category stays inside its vocabulary", f"{off} outside")

# -- did the sources actually merge ----------------------------------------

sets = db.sets()
check(len(sets) > 0, "at least one compound set was recorded")

no_set = count(
    "SELECT COUNT(*) FROM compound c WHERE NOT EXISTS ("
    "  SELECT 1 FROM compound_set_member m WHERE m.inchikey = c.inchikey)"
)
check(no_set == 0, "every compound belongs to at least one set", f"{no_set} do not")

# this is the point of merging: the same structure arriving from two sources
# has to land on one row, not two
shared = count(
    "SELECT COUNT(*) FROM (SELECT inchikey FROM compound_set_member "
    "GROUP BY inchikey HAVING COUNT(*) > 1)"
)
check(shared > 0, "compounds are shared between sets", f"{shared} in more than one set")

# and the same for targets: a protein measured by two sources is one target
multi_source = count(
    """
    SELECT COUNT(*) FROM (
      SELECT b.target_id FROM bioactivity b
        JOIN bioactivity_source s ON s.source_id = b.source_id
      GROUP BY b.target_id HAVING COUNT(DISTINCT s.source_db) > 1)
"""
)
check(
    multi_source > 0,
    "targets carry measurements from more than one source database",
    f"{multi_source} do",
)

# this is the merge test on the target side. a target is identified by the set
# of accessions it is made of, so one protein arriving from six sources has to
# end up as one row. two single-accession protein rows carrying the same
# accession would mean it did not
split = count(
    """
    SELECT COUNT(*) FROM (
      SELECT tu.uniprot_id FROM target_uniprot tu
        JOIN target t ON t.target_id = tu.target_id
       WHERE t.type = 'protein' AND t.target_id IN (
             SELECT target_id FROM target_uniprot GROUP BY target_id
              HAVING COUNT(*) = 1)
       GROUP BY tu.uniprot_id HAVING COUNT(*) > 1)
"""
)
check(split == 0, "one protein is one target however many sources report it",
      f"{split} accessions ended up on two rows")

# -- what the sources look like, which is their business and not ours -------

# a `protein` target with several accessions is a source calling a complex, an
# ortholog group or a PROTAC pair a protein. it loads, it just does not mean
# what the type says
lumped = db.read(
    """
    SELECT s.name AS set_name, COUNT(DISTINCT t.target_id) AS n
      FROM target t
      JOIN bioactivity b ON b.target_id = t.target_id
      JOIN compound_set_member m ON m.inchikey = b.inchikey
      JOIN compound_set s ON s.set_id = m.set_id
     WHERE t.type = 'protein' AND t.target_id IN (
           SELECT target_id FROM target_uniprot GROUP BY target_id HAVING COUNT(*) > 1)
     GROUP BY 1 ORDER BY n DESC
"""
)
note("protein targets carrying more than one accession, by set",
     lumped.to_string(index=False))

# a homodimer is two copies of one accession, and a set of accessions cannot
# say how many copies. so it is stored, correctly, as one member
homodimers = count(
    """
    SELECT COUNT(*) FROM (
      SELECT t.target_id FROM target t
        JOIN target_uniprot tu ON tu.target_id = t.target_id
       WHERE t.type IN ('complex', 'family', 'ppi')
       GROUP BY t.target_id HAVING COUNT(*) < 2)
"""
)
note(f"complexes with one distinct accession: {homodimers}",
     "homodimers, mostly. the schema stores a set of accessions, so it cannot "
     "record stoichiometry")

# a number with no unit and no endpoint cannot be compared with anything
blind = db.read(
    """
    SELECT s.source_db, COUNT(*) AS n,
           SUM(CASE WHEN b.bioactivity_type IS NULL THEN 1 ELSE 0 END) AS no_type
      FROM bioactivity b LEFT JOIN bioactivity_source s ON s.source_id = b.source_id
     WHERE b.unit IS NULL OR b.unit = ''
     GROUP BY 1 ORDER BY n DESC
"""
)
note(f"measurements with no unit: {int(blind['n'].sum()) if len(blind) else 0}",
     blind.to_string(index=False) if len(blind) else "")

note(f"measurements with no value: "
     f"{count('SELECT COUNT(*) FROM bioactivity WHERE value IS NULL')}",
     "the pair was tested, the number did not come through")

# a log endpoint carrying a molar unit is the one error that cannot be caught
# later: "Log Ki = 8.84 uM" reads as 8.84 uM and means 10 ** -8.84 M, which is
# 1.4 nM. six orders of magnitude, and nothing about the row looks wrong
marks = ",".join("?" * len(CONCENTRATION_UNITS))
logged = db.read(
    f"""
    SELECT bioactivity_type, unit, COUNT(*) AS n, MIN(value) AS lowest,
           MAX(value) AS highest
      FROM bioactivity
     WHERE (bioactivity_type LIKE '%log%' OR bioactivity_type LIKE 'p%50')
       AND unit IN ({marks})
     GROUP BY 1, 2 ORDER BY n DESC
""",
    *CONCENTRATION_UNITS,
)
check(
    len(logged) == 0,
    "no log endpoint carries a molar unit",
    f"{int(logged['n'].sum()) if len(logged) else 0} rows\n"
    + logged.to_string(index=False),
)

# negative concentrations would be wrong. negative % inhibition, delta Tm and
# delta G are not, so the check has to know which units it is looking at
negative = count(
    f"SELECT COUNT(*) FROM bioactivity WHERE value < 0 AND unit IN ({marks})",
    *CONCENTRATION_UNITS,
)
check(negative == 0, "no negative concentration was stored", f"{negative} rows")

note("negative values by unit, none of them a concentration",
     db.read("SELECT unit, bioactivity_type, COUNT(*) n, MIN(value) lowest "
             "  FROM bioactivity WHERE value < 0 GROUP BY 1, 2 "
             " ORDER BY n DESC LIMIT 8").to_string(index=False))

# -- report ----------------------------------------------------------------

logger.info("\n" + db.counts().to_string(index=False))
logger.info("\n" + sets.to_string(index=False))

overlap = db.read(
    """
    SELECT a.name AS set_a, b.name AS set_b, COUNT(*) AS shared
      FROM compound_set_member ma JOIN compound_set a ON a.set_id = ma.set_id
      JOIN compound_set_member mb ON mb.inchikey = ma.inchikey
      JOIN compound_set b ON b.set_id = mb.set_id
     WHERE a.set_id < b.set_id
     GROUP BY a.set_id, b.set_id ORDER BY shared DESC
"""
)
logger.info("compounds shared between sets\n" + overlap.to_string(index=False))

logger.info(
    "measurements by source database\n"
    + db.read(
        "SELECT s.source_db, COUNT(*) n FROM bioactivity b "
        "  JOIN bioactivity_source s ON s.source_id = b.source_id "
        " GROUP BY 1 ORDER BY n DESC"
    ).to_string(index=False)
)

def block(text, detail):
    body = "\n".join("        " + line for line in str(detail).splitlines())
    return text + (f"\n{body}" if detail else "")


logger.info(
    "what the sources look like, none of it a reason to distrust the merge\n"
    + "\n".join(block(f"    {what}", detail) for what, detail in notes)
)

failed = [c for c in checks if not c[0]]
logger.info(
    "\n"
    + "\n".join(
        block(f"    {'ok  ' if ok else 'FAIL'}  {what}", detail)
        for ok, what, detail in checks
    )
)
logger.info(f"{len(checks) - len(failed)}/{len(checks)} checks passed")
raise SystemExit(1 if failed else 0)
