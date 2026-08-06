"""Regressions from the adversarial review of preprocess.py.

Each test here is a defect that got past test_values.py and test_tables.py, so
each one names what it would have caught.
"""

import csv
import importlib.util
import json
import pathlib

import pytest

HERE = pathlib.Path(__file__).resolve().parent
EXPORT_JSON = HERE / "ChemicalProbesPortal-6_8_2026.json"


def load_module(path, name):
    spec = importlib.util.spec_from_file_location(name, path)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


preprocess = load_module(HERE / "preprocess.py", "preprocess")
slow = pytest.mark.slow


# ------------------------------------------------- A1: ranges and error bars

RANGES = [
    ("10-50 nM", 10.0, 50.0),
    ("0.16-0.4 uM", 0.16, 0.4),
    ("115–2000 nM", 115.0, 2000.0),
    ("5-550nM", 5.0, 550.0),
]


@pytest.mark.parametrize("raw,low,high", RANGES)
def test_a_range_is_not_written_as_a_bare_point_value(raw, low, high):
    """A range stored as its low end with no relation claims an exact potency.

    '10-50 nM' written as value=10, relation=NULL is indistinguishable from a
    measurement of exactly 10 nM. Either the operator says it is censored or the
    upper end survives somewhere, and preferably both.
    """
    fragment = preprocess.parse_potency_value(raw)[0]
    assert fragment["value"] == low
    assert fragment["value_high"] == high
    row = preprocess.bioactivity_row_from(fragment, description="ITC")
    assert row["relation"] == ">=", (
        "a range's low end has to be censored, or it reads as exact"
    )
    assert raw.replace("–", "-") in row["assay_description"].replace("–", "-"), (
        "the upper end has to survive in the text, it has no column"
    )


ERRORS = [("1 ± 0.2 nM", 1.0, 0.2), ("16 ± 11 nM", 16.0, 11.0)]


@pytest.mark.parametrize("raw,value,error", ERRORS)
def test_a_standard_deviation_survives_in_the_description(raw, value, error):
    fragment = preprocess.parse_potency_value(raw)[0]
    assert (fragment["value"], fragment["error"]) == (value, error)
    row = preprocess.bioactivity_row_from(fragment, description="SPR")
    assert "±" in row["assay_description"], (
        "the spread has no numeric column, so it lives in the text or nowhere"
    )


@slow
def test_every_range_and_error_reaches_the_written_file(tmp_path):
    """The counts behind A1: 85 range highs and 272 error bars were on no file."""
    preprocess.preprocess(EXPORT_JSON, tmp_path)
    rows = read(tmp_path, "bioactivity") + read(tmp_path, "rejected_record")
    text = " ".join(r.get("assay_description", "") + " " + r.get("label", "")
                    + " " + r.get("raw", "") for r in rows)
    probes = preprocess.load_export(EXPORT_JSON)
    missing = []
    for probe in probes:
        if not probe["InChIkey"].strip():
            continue
        for target in probe["primary_targets"]:
            for key in ("inVitroValidations", "inCellValidations"):
                for record in target.get(key) or []:
                    raw = preprocess.as_text(record.get("potencyValue"))
                    for f in preprocess.parse_potency_value(raw, record.get("potency")):
                        if f["value_high"] is None and f["error"] is None:
                            continue
                        # what has to survive is the fragment carrying the spread
                        if f["fragment"] not in text:
                            missing.append(f["fragment"])
    assert missing == [], f"{len(missing)} spreads reach no file, e.g. {missing[:5]}"


# ------------------------------------- A2: a duration is not a concentration

def test_a_time_is_never_lifted_into_the_concentration_columns():
    """'88% @ 20 h' put 20 hours into bioactivity.concentration."""
    fragment = preprocess.parse_potency_value("88% @ 20 h", "Dmax")[0]
    assert fragment["value"] == 88.0
    assert fragment["concentration_unit"] not in ("h", "min"), (
        "concentration is a concentration; a time window is not one"
    )


@pytest.mark.parametrize("raw,value,unit", [
    ("6.1 nM @ 10 uM ATP", 10.0, "uM"),
    ("35 nM (@250 uM ATP)", 250.0, "uM"),
])
def test_a_real_assay_concentration_is_still_read(raw, value, unit):
    fragment = preprocess.parse_potency_value(raw)[0]
    assert (fragment["concentration"], fragment["concentration_unit"]) == (value, unit)


# ------------------------------------------------------------- B1: the dose

def test_a_dose_does_not_borrow_the_route_of_the_dose_before_it():
    """'1 mg/Kg IV, 5 mg/Kg' gave the second dose the first one's route."""
    doses = preprocess.parse_dose("1 mg/Kg IV, 5 mg/Kg")
    assert doses[0] == (1.0, "mg/kg", "IV")
    assert doses[1][2] in (None, ""), (
        f"the 5 mg/kg dose has no stated route, got {doses[1][2]!r}"
    )


def test_a_dose_range_keeps_its_low_end_like_a_potency_does():
    """'5-30 mg/Kg PO' kept only 30, the opposite convention from bioactivity."""
    doses = preprocess.parse_dose("5-30 mg/Kg PO")
    assert [d[0] for d in doses] == [5.0], f"got {doses}"


def test_the_route_vocabulary_is_one_spelling():
    oral = {preprocess.parse_dose(f"1 mg/Kg {word}")[0][2]
            for word in ("PO", "po", "oral", "ORAL")}
    assert len(oral) == 1, f"the same route written several ways: {oral}"


# ------------------------------------------ D: report() has to actually check

def test_report_notices_a_row_that_is_both_written_and_held_back():
    """report()'s only self-check was keyed on columns a bioactivity row has
    not got, so it could never fire and main() could never exit non-zero."""
    row = {"inchikey": "A" * 14 + "-" + "B" * 10 + "-N", "target_key": "P00533",
           "value": "1", "unit": "nM", "assay_description": "SPR"}
    held = dict(row, reason="rate constant, not a potency")
    summary = preprocess.report({"bioactivity": [row], "compound": [], "target": [],
                                 "uniprot": [], "unsuitable": []},
                                {"skipped_compound": []}, [held])
    assert summary["problems"], "a row written and held back at once is a bug"


def test_report_compares_against_the_probes_it_was_given():
    summary = preprocess.report({"compound": [{"inchikey": "x"}], "bioactivity": []},
                                {"skipped_compound": []}, [], probes=[1, 2, 3])
    assert summary["problems"], "1 compound + 0 skipped is not 3 probes"


# ------------------------------- the vocabularies live in schema.sql, not here

def test_the_relations_written_come_from_the_schema():
    """README: the closed vocabularies live in the CHECK constraints in
    schema.sql and nowhere else, so they cannot drift."""
    from probedb.schema import vocabulary

    assert set(preprocess.RELATION.values()) <= vocabulary("relation")
    assert preprocess.TARGET_TYPE in vocabulary("type")


@slow
def test_no_written_relation_is_outside_the_schema_vocabulary(tmp_path):
    from probedb.schema import vocabulary

    preprocess.preprocess(EXPORT_JSON, tmp_path)
    written = {r["relation"] for r in read(tmp_path, "bioactivity") if r["relation"]}
    assert written <= vocabulary("relation")


# --------------------------------------------------- a duplicate compound key

def test_two_probes_that_resolve_to_one_key_are_not_merged_silently():
    """The resolution cache is keyed on the name, so two keyless probes with the
    same name would both take the same structure and one would vanish."""
    key = "IIIIIIIIIIIIII-JJJJJJJJJJ-N"
    probe = pd_frame([
        {"probe_ix": 0, "name": "Alpha", "inchikey": key, "smiles": "CC"},
        {"probe_ix": 1, "name": "Beta", "inchikey": key, "smiles": "CCC"},
    ])
    rows = preprocess.build_compound(probe, pd_frame([], ["probe_ix", "chembl_id"]))
    assert len({r["inchikey"] for r in rows}) == len(rows), (
        "two compound rows share an InChIKey, the loader keeps only the first"
    )


# ----------------------------------------------------------------- helpers

def read(directory, name):
    with open(pathlib.Path(directory) / f"{preprocess.PREFIX}{name}.tsv", newline="") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def pd_frame(rows, columns=None):
    import pandas as pd

    return pd.DataFrame(rows, columns=columns) if columns is not None \
        else pd.DataFrame(rows)


# ------------------------------------------- the files have to match the tables

CONTRACT = {"compound": "compound", "target": "target", "uniprot": "uniprot",
            "bioactivity": "bioactivity"}
# a column that is not a column of its own table goes somewhere the loader puts it
ELSEWHERE = {"target_key": "the join key inside the directory, never in the db",
             "chembl_id": "chembl", "source_db": "bioactivity_source",
             "source": "bioactivity_source", "xref_id": "bioactivity_source"}


@slow
@pytest.mark.parametrize("name,table", sorted(CONTRACT.items()))
def test_every_written_column_is_a_column_of_its_table(name, table, staged):
    from probedb import ProbeDB

    db = ProbeDB(":memory:", create=True)
    columns = {r[1] for r in db.conn.execute(f"PRAGMA table_info({table})")}
    db.close()
    for column in read(staged, name)[0]:
        assert column in columns or column in ELSEWHERE, (
            f"{name}.tsv writes {column!r}, which is not a column of {table} "
            f"and has nowhere else to go"
        )


@slow
def test_what_is_in_the_tables_is_what_is_in_the_files(staged):
    """Nothing invented, nothing mangled: every row in a table came from a file
    row, and every column value survives the trip."""
    db = load_staged(staged)
    blank = lambda v: "" if v is None else str(v)

    for name, sql, keys in [
        ("compound", "SELECT inchikey, smiles, name FROM compound",
         ("inchikey", "smiles", "name")),
        ("uniprot", "SELECT uniprot_id, hgnc, species FROM uniprot",
         ("uniprot_id", "hgnc", "species")),
        ("target", "SELECT type, name FROM target", ("type", "name")),
    ]:
        in_file = {tuple(r[k] for k in keys) for r in read(staged, name)}
        in_table = {tuple(blank(v) for v in row) for row in db.conn.execute(sql)}
        assert in_file == in_table, f"{name}.tsv and the table disagree"

    # chembl is unpacked from compound.tsv's '|'-joined cell
    in_file = {(c.strip().upper(), r["inchikey"]) for r in read(staged, "compound")
               for c in r["chembl_id"].split("|") if c.strip()}
    assert in_file == set(db.conn.execute("SELECT chembl_id, inchikey FROM chembl"))

    in_file = {(r["source_db"], r["source"], r["xref_id"])
               for r in read(staged, "bioactivity")}
    in_table = {tuple(blank(v) for v in row) for row in
                db.conn.execute("SELECT source_db, source, xref_id FROM bioactivity_source")}
    assert in_file == in_table, "bioactivity.tsv and bioactivity_source disagree"
    db.close()


@slow
def test_the_only_rows_missing_from_bioactivity_are_the_ones_the_loader_drops(staged):
    from collections import Counter

    db = load_staged(staged)
    blank = lambda v: "" if v is None else str(v)
    number = lambda v: float(v) if v not in (None, "") else None
    columns = ("moa", "bioactivity_type", "relation", "unit", "assay_type",
               "assay_description", "cell_line", "concentration_unit", "source_xref")
    accession = dict(db.conn.execute(
        "SELECT t.name, tu.uniprot_id FROM target_uniprot tu "
        " JOIN target t ON t.target_id = tu.target_id"))
    in_file = Counter(
        (r["inchikey"], r["target_key"], *(r[c] for c in columns),
         number(r["value"]), number(r["concentration"]))
        for r in read(staged, "bioactivity"))
    in_table = Counter(
        (row[0], accession.get(row[1], row[1]), *(blank(v) for v in row[2:11]),
         row[11], row[12])
        for row in db.conn.execute("""
            SELECT b.inchikey, t.name, b.moa, b.bioactivity_type, b.relation, b.unit,
                   b.assay_type, b.assay_description, b.cell_line,
                   b.concentration_unit, b.source_xref, b.value, b.concentration
              FROM bioactivity b JOIN target t ON t.target_id = b.target_id"""))
    assert sum((in_table - in_file).values()) == 0, (
        "the table holds rows that are in no file"
    )
    dropped = sum((in_file - in_table).values())
    assert dropped == 33, (
        f"{dropped} file rows are missing from the table, expected the 33 the "
        f"loader collapses on its 7-column identity"
    )
    db.close()


@pytest.fixture(scope="module")
def staged(tmp_path_factory):
    out = tmp_path_factory.mktemp("staged")
    preprocess.preprocess(EXPORT_JSON, out)
    return out


def load_staged(directory):
    import tempfile

    from loader import load
    from probedb import ProbeDB

    db = ProbeDB(":memory:", create=True)
    with tempfile.TemporaryDirectory() as tmp:
        load(db, preprocess.contract_directory(directory, tmp),
             source=preprocess.SOURCE_DB, strict=True)
    return db