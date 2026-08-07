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
check(split == 0,
      "no accession is stored as two separate single-accession proteins",
      f"{split} accessions ended up on two rows")

# the check above is narrow on purpose: it is the part the loader controls.
# it says nothing about an accession that also appears inside a lumped target,
# or about a target carrying no accession at all, and both of those split one
# protein across several rows just as effectively
fragmented = db.read(
    """
    SELECT tu.uniprot_id, u.hgnc, COUNT(DISTINCT tu.target_id) AS targets
      FROM target_uniprot tu
      JOIN target t ON t.target_id = tu.target_id
      JOIN uniprot u ON u.uniprot_id = tu.uniprot_id
     WHERE t.type = 'protein'
     GROUP BY tu.uniprot_id HAVING targets > 1
     ORDER BY targets DESC LIMIT 6
"""
)
total_fragmented = count(
    """
    SELECT COUNT(*) FROM (
      SELECT tu.uniprot_id FROM target_uniprot tu
        JOIN target t ON t.target_id = tu.target_id
       WHERE t.type = 'protein'
       GROUP BY tu.uniprot_id HAVING COUNT(DISTINCT tu.target_id) > 1)
"""
)
note(f"accessions reached by more than one protein target: {total_fragmented}",
     fragmented.to_string(index=False))

# db.targets_for() matches an HGNC symbol as well as an accession, so a symbol
# claimed by two accessions silently answers for the wrong protein. P01112 is
# HRAS and arrives here filed under hgnc 'KRAS', so asking for KRAS returns the
# HRAS target too. an accession is unambiguous, a symbol is not
ambiguous = db.read(
    """
    SELECT hgnc, COUNT(*) AS accessions, GROUP_CONCAT(uniprot_id) AS which
      FROM uniprot WHERE hgnc IS NOT NULL AND hgnc != ''
     GROUP BY hgnc HAVING accessions > 1
     ORDER BY accessions DESC LIMIT 6
"""
)
total_ambiguous = count(
    """
    SELECT COUNT(*) FROM (
      SELECT hgnc FROM uniprot WHERE hgnc IS NOT NULL AND hgnc != ''
       GROUP BY hgnc HAVING COUNT(*) > 1)
"""
)
note(f"HGNC symbols claimed by more than one accession: {total_ambiguous}",
     ambiguous.to_string(index=False)
     + "\nmost are orthologs sharing a gene name, but KRAS is not: it is on "
       "P01116 and on P01112, which is HRAS. select on the accession when it "
       "matters")

# a target with no accession is identified by its name alone, so it can only
# merge with a source that spells the name identically. EGFR arrives as
# 'EGFR' with an accession and as 'Epidermal growth factor receptor' without
# one, and those are two targets that will never come together
nameless = count(
    """
    SELECT COUNT(*) FROM target t WHERE NOT EXISTS
      (SELECT 1 FROM target_uniprot tu WHERE tu.target_id = t.target_id)
"""
)
nameless_rows = count(
    """
    SELECT COUNT(*) FROM bioactivity b WHERE NOT EXISTS
      (SELECT 1 FROM target_uniprot tu WHERE tu.target_id = b.target_id)
"""
)
note(f"targets with no accession at all: {nameless}",
     f"they carry {nameless_rows} measurements and can only ever merge with a "
     f"source that spells the name the same way")

# -- what the sources look like, which is their business and not ours -------

# a `protein` target with several accessions is a source calling a complex, an
# ortholog group or a whole screening panel a protein. it loads, it just does
# not mean what the type says, and now that uniprot.superfamily exists it does
# real damage: every one of those accessions drags the target into its family,
# so asking for the Ras family also returns whatever else got lumped in
lumped = db.read(
    """
    SELECT t.target_id, t.name, COUNT(*) AS accessions
      FROM target t JOIN target_uniprot tu ON tu.target_id = t.target_id
     WHERE t.type = 'protein'
     GROUP BY t.target_id HAVING accessions > 1
     ORDER BY accessions DESC LIMIT 8
"""
)
total = count(
    """
    SELECT COUNT(*) FROM (
      SELECT tu.target_id FROM target_uniprot tu
        JOIN target t ON t.target_id = tu.target_id
       WHERE t.type = 'protein'
       GROUP BY tu.target_id HAVING COUNT(*) > 1)
"""
)
note(f"protein targets carrying more than one accession: {total}",
     lumped.to_string(index=False))

# how much of the reference file landed, and how far it reaches
classified = count("SELECT COUNT(*) FROM uniprot WHERE superfamily IS NOT NULL")
reached = count(
    """
    SELECT COUNT(*) FROM bioactivity b WHERE b.target_id IN (
      SELECT tu.target_id FROM target_uniprot tu
        JOIN uniprot u ON u.uniprot_id = tu.uniprot_id
       WHERE u.superfamily IS NOT NULL)
"""
)
measurements = count("SELECT COUNT(*) FROM bioactivity")
note(f"accessions with a UniProt family: {classified} of "
     f"{count('SELECT COUNT(*) FROM uniprot')}",
     f"they carry {reached} of {measurements} measurements "
     f"({100 * reached // max(measurements, 1)}%), so the classified proteins "
     f"are the well studied ones")

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
