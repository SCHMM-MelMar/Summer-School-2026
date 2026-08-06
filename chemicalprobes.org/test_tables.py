"""Tests for preprocess.py: the table builders, the writer, and the load.

Written against the contract, which is the docstring of each function plus
README section 3 and the counts in README.md -- not against a reading of
the bodies. A test that fails here is either a function that does not keep its
docstring or a count that has moved.

What is covered: load_export, flatten, build_compound, build_target,
build_uniprot, build_leftovers,
resolve_missing_structures, write_table, write_staging, report and main. The
value parsers (parse_potency_value, canonical_unit, repair, ...) have their own
suite; they appear here only where a table builder has to carry their result.

The tests that read the real export are marked `slow` and named
`test_integration_*`, so either of

    pytest test_tables.py -m "not slow"
    pytest test_tables.py -k "not integration"

leaves the unit tests. The marker is not registered anywhere: this suite adds no
pytest.ini and no conftest.py by instruction, so pytest warns once that `slow`
is unknown. The warning is harmless.

Nothing here touches the network, and nothing here writes into the repository.
"""

import csv
import importlib.util
import inspect
import json
import math
import pathlib
import re
import socket
from collections.abc import Mapping

import pandas as pd
import pytest
from rdkit import Chem
from rdkit.Chem.inchi import MolToInchiKey

from loader import load as _load
from loader import validate
from loader.validate import read as _read
from probedb import ProbeDB
from probedb.schema import vocabulary

slow = pytest.mark.slow

# what "fails loudly" means for load_export. NotImplementedError is a
# RuntimeError and deliberately not in this tuple, so a function that is only a
# stub cannot pass a test that asks it to reject something.
SHAPE_ERRORS = (ValueError, KeyError, TypeError, AssertionError)

HERE = pathlib.Path(__file__).resolve().parent
PREPROCESS_PY = HERE / "preprocess.py"
EXPORT_JSON = HERE / "ChemicalProbesPortal-6_8_2026.json"

# the cache of externally resolved structures that sits next to preprocess.py,
# with the columns it has today. D1 says the lookup is done once and recorded;
# the tests read that file and never a network.
RESOLVED_CACHE = "resolved_structures.tsv"
CACHE_COLUMNS = ["name", "inchikey", "smiles", "resolved_from", "lookup",
                 "resolved_title"]


def load_module(path, name):
    spec = importlib.util.spec_from_file_location(name, path)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


# the directory name contains a dot, so preprocess is not importable as a package
preprocess = load_module(PREPROCESS_PY, "preprocess")


# ------------------------------------------------------------- the contract

# README section 3, and the docstrings of the builders
COMPOUND_COLUMNS = ["inchikey", "smiles", "chembl_id", "name"]
TARGET_COLUMNS = ["target_key", "type", "name"]
UNIPROT_COLUMNS = ["uniprot_id", "target_key", "hgnc", "species"]
BIOACTIVITY_COLUMNS = [
    "inchikey", "target_key", "moa", "bioactivity_type", "relation", "value",
    "unit", "assay_type", "assay_description", "cell_line", "concentration",
    "concentration_unit", "source_db", "source", "source_xref", "xref_id",
]
BIOACTIVITY_REQUIRED = ["inchikey", "target_key", "value", "unit"]
ASSESSMENT_COLUMNS = [
    "inchikey", "verdict", "pains", "toxicophore", "rating_in_cell",
    "rating_in_organism", "rating_count", "published_date", "source_db", "source",
]

# the 18 keys of a probe record
EXPORT_KEYS = [
    "name", "rating_in_cell", "rating_in_organism", "rating_count",
    "unsuitable", "URL", "primary_targets", "published_date", "InChIkey",
    "smiles", "pains", "toxicophore", "canSAR_ID", "ChEMBL_ID", "PMID",
    "inVitroValidations", "inVivoValidations", "control_compounds",
]

FRAME_NAMES = ("probe", "target", "validation", "invivo", "chembl",
               "reference", "control")

# README.md, "What is in the export"
EXPORT_FRAME_ROWS = {
    "probe": 1247,
    "target": 1372,
    "validation": 3551,
    "invivo": 1265,
    "chembl": 678,
    "reference": 1816,
    "control": 418,
}

# README.md, "What would be written"
EXPECTED_COMPOUNDS = 1223
EXPECTED_TARGETS = 644
EXPECTED_UNIPROTS = 644
EXPECTED_BIOACTIVITY_NUMBERS = 3649
EXPECTED_QUARANTINED = 23
EXPECTED_UNSUITABLE = 260
EXPECTED_SKIPPED_COMPOUNDS = 24
EXPECTED_PROBES = 1247
# the loader is not touched by decision, so this drop is current behaviour
EXPECTED_DUPLICATES_DROPPED = 33   # 35 before ranges were censored with '>='
EXPECTED_GENUINE_REPEATS = 6
# 156 in README.md, 149 measured: see the pinned reload test
EXPECTED_REINSERTED_ON_RELOAD = 149


# ---------------------------------------------------------------- fixtures

NOTE = "Please cite the Chemical Probes Portal"

# assayDesc carries tabs, line breaks and double quotes. No leading or trailing
# whitespace, because loader.validate.read() strips every value it reads back.
HAZARD_DESC = (
    'SPR direct binding assay\trun in duplicate,\n"cell-free" throughout'
)

PROBE_WITH_NO_POTENCY_KEY = {
    "name": "SGC0946",
    "rating_in_cell": 3.75,
    "rating_in_organism": 0,
    "rating_count": 4,
    "unsuitable": "No",
    "URL": "https://www.chemicalprobes.org/sgc0946",
    "primary_targets": [
        {
            "name": "DOT1L",
            "class": "Epigenetic",
            "subClass": "Protein methyltransferase",
            "moa": "Inhibitor",
            # an in-vitro validation has no `potency` key at all
            "inVitroValidations": [
                {"potencyValue": "0.06 nM", "assayDesc": HAZARD_DESC}
            ],
            "inCellValidations": [],
            "uniprot_id": "Q8TEK3",
        }
    ],
    "published_date": "2015-10-02",
    "InChIkey": "IQCKJUKAQJINMK-HUBRGWSESA-N",
    "smiles": "CC(C)N(CCCNC(=O)Nc1ccc(C(C)(C)C)cc1)C[C@H]1O[C@@H](n2cc(Br)c3c"
              "(N)ncnc32)[C@H](O)[C@@H]1O",
    "pains": "No",
    "toxicophore": "No",
    "canSAR_ID": 828356,
    "ChEMBL_ID": ["CHEMBL3087498"],
    "PMID": [
        "www.doi.org/10.1038/ncomms2304",
        "http://www.ncbi.nlm.nih.gov/pubmed/23250418",
    ],
    "inVitroValidations": [],
    "inVivoValidations": [],
    "control_compounds": ["SGC0649"],
}

PROBE_WITH_TWO_CHEMBL_IDS = {
    "name": "MZ1 ",                       # 10 names carry whitespace
    "rating_in_cell": 3.0,
    "rating_in_organism": 0,
    "rating_count": 2,
    "unsuitable": "No",
    "URL": "https://www.chemicalprobes.org/mz1",
    "primary_targets": [
        {
            "name": "BRD4",
            "class": "Epigenetic",
            "subClass": "Bromodomain",
            "moa": "Degrader (PROTAC)",
            "inVitroValidations": [],
            # one potencyValue holding two measurements, each naming its endpoint
            "inCellValidations": [
                {
                    "potency": "DC50, Dmax",
                    "potencyValue": "63 nM (DC50); 90.8% (Dmax)",
                    "assayDesc": "degradation of BRD4 in MV4;11 cells",
                }
            ],
            "uniprot_id": "O60885",
        }
    ],
    "published_date": "2017-05-11",
    "InChIkey": "XZXHXSATPCNXJR-ZIADKAODSA-N",
    "smiles": "CC(C)Nc1ncnc2c1ncn2C",
    "pains": "No",
    "toxicophore": "No",
    "canSAR_ID": 1013442,
    "ChEMBL_ID": ["CHEMBL3545209", "CHEMBL4297397"],
    "PMID": [],
    "inVitroValidations": [],
    "inVivoValidations": [],
    "control_compounds": [],
}

PROBE_UNSUITABLE = {
    "name": "JIB-04",
    "rating_in_cell": 0,
    "rating_in_organism": 0,
    "rating_count": 0,
    "unsuitable": "Yes",
    # 260 unsuitable probes live under /unsuitables/
    "URL": "https://www.chemicalprobes.org/unsuitables/jib-04",
    "primary_targets": [],                # no target, no validation, no in vivo
    "published_date": "2016-12-15",
    "InChIkey": "YHHFKWKMXWRVTJ-OQKWZONESA-N",
    "smiles": "Clc1ccc(N/N=C(\\c2ccccc2)c2ccccn2)nc1",
    "pains": "No",
    "toxicophore": "Yes",
    "canSAR_ID": 1354531,
    "ChEMBL_ID": [],
    "PMID": ["www.doi.org/10.1000/example"],
    "inVitroValidations": [],
    "inVivoValidations": [],
    "control_compounds": [],
}

PROBE_WITH_NO_INCHIKEY = {
    "name": "CBK-1234",                   # an internal code, not in the cache
    "rating_in_cell": 2.0,
    "rating_in_organism": 0,
    "rating_count": 1,
    "unsuitable": "No",
    "URL": "https://www.chemicalprobes.org/cbk-1234",
    "primary_targets": [
        {
            "name": "EGFR",
            "class": "Kinase",
            "subClass": "Protein kinase",
            "moa": "Inhibitor",
            "inVitroValidations": [],
            "inCellValidations": [
                {
                    "potency": "IC50",
                    "potencyValue": "12 nM",
                    "assayDesc": "phospho-EGFR in A431 cells",
                }
            ],
            "uniprot_id": "P00533",
        }
    ],
    "published_date": "2024-02-01",
    "InChIkey": "",
    "smiles": "",
    "pains": "No",
    "toxicophore": "No",
    "canSAR_ID": 999001,
    "ChEMBL_ID": [],
    "PMID": [],
    "inVitroValidations": [],
    "inVivoValidations": [],
    "control_compounds": [],
}

PROBE_WITH_IN_VIVO = {
    "name": "AT7519",
    "rating_in_cell": 3.5,
    "rating_in_organism": 3.0,
    "rating_count": 3,
    "unsuitable": "No",
    "URL": "https://www.chemicalprobes.org/at7519",
    "primary_targets": [
        {
            "name": "RPS6KA3",
            "class": "Kinase",
            "subClass": "Protein kinase",
            "moa": "Inhibitor",
            # a reciprocal rate constant: quarantined, never a potency
            "inVitroValidations": [
                {
                    "potencyValue": "9900 ± 1800 M-1 s-1",
                    "assayDesc": "stopped-flow kinetics of adduct formation",
                }
            ],
            "inCellValidations": [],
            "uniprot_id": "P51812",
        }
    ],
    "published_date": "2018-06-30",
    "InChIkey": "AQGNHMYCWQNFFB-UHFFFAOYSA-N",
    "smiles": "CCOc1ccccc1",
    "pains": "No",
    "toxicophore": "No",
    "canSAR_ID": 705632,
    "ChEMBL_ID": [],
    "PMID": [],
    "inVitroValidations": [],
    "inVivoValidations": [
        {"organism": "Mouse", "dose": "1 mg/Kg IV, 5 mg/Kg PO"},
        {"organism": "Rat"},              # 7 in vivo records have no `dose` key
    ],
    "control_compounds": [],
}

PROBES = [
    PROBE_WITH_NO_POTENCY_KEY,
    PROBE_WITH_TWO_CHEMBL_IDS,
    PROBE_UNSUITABLE,
    PROBE_WITH_NO_INCHIKEY,
    PROBE_WITH_IN_VIVO,
]

KEYED_PROBES = [p for p in PROBES if p["InChIkey"]]

ETHANOL = "CCO"
ETHANOL_KEY = MolToInchiKey(Chem.MolFromSmiles(ETHANOL))
WRONG_KEY = "AAAAAAAAAAAAAA-BBBBBBBBBB-N"


# ----------------------------------------------------------------- helpers

def rows_of(frame):
    """A builder's return value as a list of dicts."""
    if frame is None:
        return []
    if hasattr(frame, "to_dict"):
        return frame.to_dict("records")
    if isinstance(frame, Mapping):
        return [frame]
    return list(frame)


def columns_of(frame):
    if hasattr(frame, "columns"):
        return list(frame.columns)
    rows = rows_of(frame)
    return list(rows[0]) if rows else []


def text(value):
    """What a value looks like once written: NaN and None are nothing."""
    if value is None:
        return ""
    if isinstance(value, float) and math.isnan(value):
        return ""
    return str(value)


def cell(row, column):
    return text(row.get(column))


def read_staged(directory, name):
    """loader/'s reader, on the prefixed names preprocess.py writes."""
    import tempfile

    directory = pathlib.Path(directory)
    prefixed = directory / f"{preprocess.PREFIX}{name}.tsv"
    if not prefixed.exists():
        return _read(directory, name)
    with tempfile.TemporaryDirectory() as tmp:
        (pathlib.Path(tmp) / f"{name}.tsv").write_bytes(prefixed.read_bytes())
        return _read(tmp, name)


def validate_staged(directory):
    """loader.validate on the contract names, which is what it reads."""
    import tempfile

    with tempfile.TemporaryDirectory() as tmp:
        return validate(preprocess.contract_directory(directory, tmp))


def load_staging(db, directory, **kw):
    """loader/ wants the four contract names, which is what production hands it."""
    import tempfile

    with tempfile.TemporaryDirectory() as tmp:
        return _load(db, preprocess.contract_directory(directory, tmp),
                     source=preprocess.SOURCE_DB, **kw)


def header_of(path):
    with open(path, newline="") as handle:
        return next(csv.reader(handle, delimiter="\t"))


class Flat:
    """flatten()'s seven frames, however it hands them back."""

    def __init__(self, result):
        if isinstance(result, Mapping):
            frames = dict(result)
        else:
            names = getattr(result, "_fields", FRAME_NAMES)
            values = list(result)
            assert len(values) == len(FRAME_NAMES), (
                "flatten() returns the seven frames of its docstring"
            )
            frames = dict(zip(names, values))
        missing = sorted(set(FRAME_NAMES) - set(frames))
        assert not missing, f"flatten() gives no {missing} frame"
        self.frames = frames

    def raw(self, name):
        return self.frames[name]

    def rows(self, name):
        return rows_of(self.frames[name])


def flat():
    """The fixture probes through flatten()."""
    return Flat(preprocess.flatten(PROBES))


def call(function, **available):
    """Call a builder with the arguments it declares, by name.

    Every builder names its arguments after the frame it reads, so a test offers
    all of them and lets the signature choose. A builder that also takes the
    reference frame, or the set of keyed probes, is then called the way the
    pipeline calls it, and a test does not have to track argument order.
    """
    wanted = inspect.signature(function).parameters
    return function(**{name: available[name] for name in wanted
                       if name in available})


def offer(f):
    """Every argument the builders take, under the name they use for it."""
    return {
        "probes": PROBES,
        "frames": f.frames,
        "probe": f.raw("probe"),
        "target": f.raw("target"),
        "validation": f.raw("validation"),
        "invivo": f.raw("invivo"),
        "chembl": f.raw("chembl"),
        "reference": f.raw("reference"),
        "control": f.raw("control"),
        "keyed": {r["probe_ix"] for r in f.rows("probe")
                  if text(r.get("inchikey"))},
    }


def written_and_held(arguments):
    """The bioactivity rows written, and the rows held back (D15).

    build_bioactivity() gives the rows its docstring describes. The held-back
    ones have to go somewhere too, so they are taken from a second return value
    or from the function that builds quarantine.tsv, and from the parser only if
    the module names neither.
    """
    result = call(preprocess.split_bioactivity, **arguments)
    if isinstance(result, tuple) and len(result) == 2:
        return rows_of(result[0]), rows_of(result[1])
    held = getattr(preprocess, "build_quarantine", None)
    return (rows_of(result),
            rows_of(call(held, **arguments)) if held else quarantined_numbers())


def built(f=None):
    """Every frame the writer is handed, straight from the builders."""
    f = f or flat()
    arguments = offer(f)
    rows, quarantined = written_and_held(arguments)
    tables = {
        "compound": call(preprocess.build_compound, **arguments),
        "target": call(preprocess.build_target, **arguments),
        "uniprot": call(preprocess.build_uniprot, **arguments),
        "bioactivity": rows,
    }
    # the held-back rows are one stage of rejected_record now, so they go in
    # through build_leftovers rather than being a table of their own
    leftovers = call(preprocess.build_leftovers, quarantined=quarantined, **arguments)
    return tables, leftovers


def parsed_numbers(probe):
    """Every number the parser reads out of one probe's validations."""
    out = []
    for target in probe["primary_targets"]:
        validations = (target["inVitroValidations"]
                       + target["inCellValidations"])
        for record in validations:
            out += rows_of(
                preprocess.parse_potency_value(
                    record["potencyValue"], record.get("potency")
                )
            )
    return out


def candidate_numbers():
    """Numbers under a compound that can be written at all."""
    return [n for p in KEYED_PROBES for n in parsed_numbers(p)]


def quarantined_numbers():
    return [n for n in candidate_numbers() if n.get("quarantine")]


def frame_in(result, name):
    """One named frame out of whatever build_leftovers hands back."""
    assert isinstance(result, Mapping), (
        "build_leftovers() returns its frames keyed by kind"
    )
    assert name in result, f"build_leftovers() gives no {name!r} frame"
    return rows_of(result[name])


def as_count(value):
    if isinstance(value, bool):
        return None
    if isinstance(value, int):
        return value
    if isinstance(value, Mapping):
        for key in ("rows", "count", "n", "written"):
            if isinstance(value.get(key), int):
                return value[key]
        return None
    if isinstance(value, (list, tuple)):
        return len(value)
    return None


def count_in(report, name):
    """The row count report() gives for one file, however it nests it."""
    if isinstance(report, Mapping):
        if name in report:
            return as_count(report[name])
        for value in report.values():
            found = count_in(value, name)
            if found is not None:
                return found
    return None


def loader_identity(row):
    """The duplicate key loader/load.py:75-83 builds, on a staged row.

    Mirrored here so a test can count what the loader will drop without
    changing anything in loader/.
    """
    return (
        row["inchikey"].upper(),
        row["target_key"],
        (row.get("source_db") or None, row.get("source") or None),
        row.get("bioactivity_type") or None,
        row.get("relation") or None,
        None if not row.get("value") else float(row["value"]),
        row.get("unit") or None,
    )


def block_network(monkeypatch):
    """No test may reach the network."""
    def refuse(*args, **kwargs):
        raise RuntimeError("a test tried to open a network connection")

    monkeypatch.setattr(socket.socket, "connect", refuse, raising=False)
    monkeypatch.setattr(socket.socket, "connect_ex", refuse, raising=False)
    monkeypatch.setattr(socket, "create_connection", refuse, raising=False)


@pytest.fixture
def no_network(monkeypatch):
    block_network(monkeypatch)


def module_with_cache(directory, cache_rows):
    """A copy of preprocess.py in `directory`, with the cache it reads.

    The cache is looked up next to the module, so copying the module puts the
    test in charge of it without writing anything into the repository.
    cache_rows None means no cache file at all, [] means a header and no rows.
    """
    directory.mkdir(parents=True, exist_ok=True)
    copy = directory / "preprocess.py"
    copy.write_bytes(PREPROCESS_PY.read_bytes())
    if cache_rows is not None:
        with open(directory / RESOLVED_CACHE, "w", newline="") as handle:
            writer = csv.DictWriter(handle, fieldnames=CACHE_COLUMNS,
                                    delimiter="\t")
            writer.writeheader()
            for row in cache_rows:
                writer.writerow(row)
    return load_module(copy, "preprocess_" + directory.name)


def resolve(module, probes):
    """What resolve_missing_structures() recovers, as rows.

    Its docstring names the argument `probe`, which is the flattened frame the
    other builders take, so that is what is offered.
    """
    frame = Flat(module.flatten(probes)).raw("probe")
    return rows_of(call(module.resolve_missing_structures,
                        probe=frame, probes=probes))


_BUILT = {}


@pytest.fixture(scope="session")
def workdir(tmp_path_factory):
    return tmp_path_factory.mktemp("preprocess")


@pytest.fixture(scope="session")
def fixture_export(workdir):
    path = workdir / "fixture.json"
    path.write_text(json.dumps({"Note": NOTE, "probes": PROBES}),
                    encoding="utf-8")
    return path


def build_staging(workdir, json_path, name):
    if name not in _BUILT:
        out = workdir / name
        out.mkdir(exist_ok=True)
        with pytest.MonkeyPatch.context() as patch:
            block_network(patch)
            preprocess.main(["--json", str(json_path), "--out", str(out)])
        _BUILT[name] = out
    return _BUILT[name]


# The two staging fixtures hand back a callable rather than a path, so a failure
# inside main() is reported against the test that needed the directory rather
# than as a fixture error. Each directory is built once per session.

@pytest.fixture(scope="session")
def fixture_staging(workdir, fixture_export):
    return lambda: build_staging(workdir, fixture_export, "fixture")


@pytest.fixture(scope="session")
def portal_staging(workdir):
    return lambda: build_staging(workdir, EXPORT_JSON, "portal")


@pytest.fixture(scope="session")
def portal_db(portal_staging):
    def build():
        if "db" not in _BUILT:
            db = ProbeDB(":memory:", create=True)
            report = load_staging(db, portal_staging(), strict=True)
            _BUILT["db"] = (db, report)
        return _BUILT["db"]

    return build


def export_probes():
    """The real export, read here only as reference data for assertions."""
    if "export" not in _BUILT:
        _BUILT["export"] = json.loads(EXPORT_JSON.read_text())["probes"]
    return _BUILT["export"]


def export_frames():
    if "flat" not in _BUILT:
        _BUILT["flat"] = Flat(
            preprocess.flatten(preprocess.load_export(EXPORT_JSON))
        )
    return _BUILT["flat"]


# ---------------------------------------------------------------- load_export

def test_load_export_rejects_an_unknown_key(tmp_path):
    # a new key in a later release has to stop the run, not be ignored
    payload = {"Note": NOTE, "probes": [dict(PROBES[0], newField="x")]}
    path = tmp_path / "extra.json"
    path.write_text(json.dumps(payload))
    with pytest.raises(SHAPE_ERRORS):
        preprocess.load_export(path)


@pytest.mark.parametrize("payload", [
    {},                                        # no probes at all
    {"Note": NOTE, "probes": {}},              # probes is not a list
    {"Note": NOTE, "probes": [{"name": "x"}]},  # a probe missing 17 keys
    {"Note": NOTE, "probes": [[]]},            # a probe that is not an object
])
def test_load_export_rejects_a_bad_shape(tmp_path, payload):
    path = tmp_path / "bad.json"
    path.write_text(json.dumps(payload))
    with pytest.raises(SHAPE_ERRORS):
        preprocess.load_export(path)


@slow
def test_integration_load_export_reads_every_probe():
    probes = preprocess.load_export(EXPORT_JSON)
    assert len(probes) == EXPECTED_PROBES
    assert all(isinstance(p, dict) for p in probes)
    assert all(set(p) == set(EXPORT_KEYS) for p in probes)


# -------------------------------------------------------------------- flatten

@pytest.mark.parametrize("name,rows", [
    ("probe", 5),
    ("target", 4),        # the unsuitable probe lists none
    ("validation", 4),
    ("invivo", 2),
    ("chembl", 3),        # one probe carries two ids
    ("reference", 3),
    ("control", 1),
])
def test_flatten_frame_sizes(name, rows):
    assert len(flat().rows(name)) == rows


def test_flatten_keeps_a_validation_with_no_potency_key():
    # an in-vitro validation has no `potency` key at all, so .get() throughout.
    # flatten's columns are internal to the module, so the record is found by
    # its value and not by a column name.
    rows = flat().rows("validation")
    kept = [r for r in rows if "0.06 nM" in {text(v) for v in r.values()}]
    assert len(kept) == 1
    assert "IC50" not in {text(v) for v in kept[0].values()}


def test_flatten_keeps_an_in_vivo_record_with_no_dose_key():
    rows = flat().rows("invivo")
    organisms = {cell(r, "organism") for r in rows}
    assert organisms == {"Mouse", "Rat"}
    rat = [r for r in rows if cell(r, "organism") == "Rat"][0]
    assert cell(rat, "dose") == ""


@slow
@pytest.mark.parametrize("name,rows", sorted(EXPORT_FRAME_ROWS.items()))
def test_integration_flatten_frame_sizes(name, rows):
    assert len(export_frames().rows(name)) == rows


# ------------------------------------------------------------- build_compound

def compound_rows(f=None):
    return rows_of(call(preprocess.build_compound, **offer(f or flat())))


def test_build_compound_columns_are_the_contract():
    assert columns_of(
        call(preprocess.build_compound, **offer(flat()))
    ) == COMPOUND_COLUMNS


def test_build_compound_joins_several_chembl_ids_with_a_pipe():
    # loader/load.py:52 splits chembl_id on '|'
    rows = compound_rows()
    row = [r for r in rows
           if cell(r, "inchikey") == PROBE_WITH_TWO_CHEMBL_IDS["InChIkey"]][0]
    assert cell(row, "chembl_id") == "CHEMBL3545209|CHEMBL4297397"
    assert set(cell(row, "chembl_id").split("|")) == {
        "CHEMBL3545209", "CHEMBL4297397"
    }


def test_build_compound_strips_the_name():
    names = [cell(r, "name") for r in compound_rows()]
    assert "MZ1" in names
    assert all(name == name.strip() for name in names)


def test_build_compound_leaves_out_a_probe_with_no_inchikey():
    rows = compound_rows()
    assert len(rows) == len(KEYED_PROBES)
    assert all(cell(r, "inchikey") for r in rows)
    assert PROBE_WITH_NO_INCHIKEY["name"] not in [cell(r, "name") for r in rows]


# --------------------------------------------------- build_target / uniprot

def target_rows():
    return rows_of(call(preprocess.build_target, **offer(flat())))


def uniprot_rows():
    return rows_of(call(preprocess.build_uniprot, **offer(flat())))


def test_build_target_columns_are_the_contract():
    assert columns_of(
        call(preprocess.build_target, **offer(flat()))
    ) == TARGET_COLUMNS


def test_build_target_is_one_protein_row_per_accession():
    rows = target_rows()
    keys = [cell(r, "target_key") for r in rows]
    assert len(keys) == len(set(keys))
    assert {"Q8TEK3", "O60885", "P51812"} <= set(keys)
    # class and subClass are a family taxonomy, never target.type
    assert {cell(r, "type") for r in rows} == {"protein"}
    names = {cell(r, "target_key"): cell(r, "name") for r in rows}
    assert names["Q8TEK3"] == "DOT1L"


def test_build_uniprot_columns_are_the_contract():
    assert columns_of(
        call(preprocess.build_uniprot, **offer(flat()))
    ) == UNIPROT_COLUMNS


def test_build_uniprot_keys_on_the_accession_and_leaves_species_empty():
    rows = uniprot_rows()
    assert rows
    for row in rows:
        assert cell(row, "uniprot_id") == cell(row, "target_key")
        assert cell(row, "species") == ""      # a UniProt lookup, not in the export
    hgnc = {cell(r, "uniprot_id"): cell(r, "hgnc") for r in rows}
    assert hgnc["O60885"] == "BRD4"


# --------------------------------------------------------- build_bioactivity

def bioactivity_rows():
    rows, _ = written_and_held(offer(flat()))
    return rows


def test_build_bioactivity_columns_are_the_contract():
    assert columns_of(bioactivity_rows()) == BIOACTIVITY_COLUMNS


@pytest.mark.parametrize("target_key,assay_type", [
    ("Q8TEK3", "binding"),   # in vitro, and assayDesc names SPR (D7)
    ("O60885", "cell"),      # in cell
])
def test_build_bioactivity_assay_type_per_tier(target_key, assay_type):
    rows = [r for r in bioactivity_rows()
            if cell(r, "target_key") == target_key]
    assert rows
    assert {cell(r, "assay_type") for r in rows} == {assay_type}


def test_build_bioactivity_leaves_cell_line_empty():
    # D10: a wrong line is worse than none
    assert all(cell(r, "cell_line") == "" for r in bioactivity_rows())


def test_build_bioactivity_records_the_portal_as_the_source():
    url = {p["InChIkey"]: p["URL"] for p in KEYED_PROBES}
    for row in bioactivity_rows():
        assert cell(row, "source_db") == preprocess.SOURCE_DB
        assert cell(row, "source")
        assert cell(row, "xref_id") == preprocess.PORTAL_PREFIX
        # the full path, not the last segment: /unsuitables/ would 404
        assert (cell(row, "xref_id") + cell(row, "source_xref")
                == url[cell(row, "inchikey")])


def test_build_bioactivity_takes_the_endpoint_per_fragment():
    # D9: two measurements in one potencyValue, each naming its endpoint. No row
    # may carry the comma-joined label of all of them.
    rows = [r for r in bioactivity_rows() if cell(r, "target_key") == "O60885"]
    assert len(rows) == 2
    assert {cell(r, "bioactivity_type") for r in rows} == {"DC50", "Dmax"}
    assert all("," not in cell(r, "bioactivity_type") for r in rows)
    # the fragment label travels with the number, so the two differ
    descriptions = [cell(r, "assay_description") for r in rows]
    assert descriptions[0] != descriptions[1]
    raw = PROBE_WITH_TWO_CHEMBL_IDS["primary_targets"][0]
    raw = raw["inCellValidations"][0]["assayDesc"]
    assert all(raw in d for d in descriptions)


def test_build_bioactivity_moa_is_normalised():
    # moa is part of the bioactivity_group key and the composite FK is byte-exact
    rows = [r for r in bioactivity_rows() if cell(r, "target_key") == "O60885"]
    assert {cell(r, "moa") for r in rows} == {"degrader (protac)"}


def test_build_bioactivity_writes_nothing_for_a_keyless_probe():
    rows = bioactivity_rows()
    assert "P00533" not in {cell(r, "target_key") for r in rows}
    assert all(cell(r, "inchikey") for r in rows)


# --------------------------------------------------------- build_leftovers

def test_build_leftovers_keeps_the_in_vivo_record_with_no_dose():
    leftovers = call(preprocess.build_leftovers, **offer(flat()))
    rows = frame_in(leftovers, "in_vivo_dose")
    assert {cell(r, "organism") for r in rows} == {"Mouse", "Rat"}
    rat = [r for r in rows if cell(r, "organism") == "Rat"]
    assert len(rat) == 1
    assert cell(rat[0], "dose_value") == ""
    assert all(cell(r, "inchikey") == PROBE_WITH_IN_VIVO["InChIkey"]
               for r in rows)


def test_build_leftovers_keeps_class_per_target_key():
    leftovers = call(preprocess.build_leftovers, **offer(flat()))
    rows = frame_in(leftovers, "target_class")
    dot1l = [r for r in rows if cell(r, "target_key") == "Q8TEK3"]
    assert dot1l
    values = {text(v) for row in dot1l for v in row.values()}
    assert {"Epigenetic", "Protein methyltransferase"} <= values


def test_build_leftovers_keeps_the_controls_and_the_url_per_compound():
    leftovers = call(preprocess.build_leftovers, **offer(flat()))
    rows = frame_in(leftovers, "compound_annotation")
    key = PROBE_WITH_NO_POTENCY_KEY["InChIkey"]
    mine = [r for r in rows if cell(r, "inchikey") == key]
    values = {text(v) for row in mine for v in row.values()}
    assert "SGC0649" in values                            # a control compound
    # the url is not an annotation any more: it is the source, and the flags and
    # dates moved to probe_assessment, so this table is the tail it should be
    assert {text(r.get("property")) for r in rows} <= {"control_compound",
                                                      "structure_source"}
    written = {p["InChIkey"] for p in KEYED_PROBES}
    assert {cell(r, "inchikey") for r in rows} <= written


# ------------------------------------------------ resolve_missing_structures

@pytest.mark.parametrize("cache", ["absent", "header-only", "zero-bytes"])
def test_resolve_returns_nothing_without_a_cache(tmp_path, no_network, cache):
    directory = tmp_path / cache.replace("-", "_")
    module = module_with_cache(directory, None if cache == "absent" else [])
    if cache == "zero-bytes":
        (directory / RESOLVED_CACHE).write_text("")
    assert resolve(module, PROBES) == []


def test_resolve_does_not_reach_the_network(tmp_path, no_network):
    # the lookup was done once and recorded; a test run resolves from the cache
    # or not at all
    module = module_with_cache(tmp_path / "no_row", [
        {"name": "SOMETHING ELSE", "inchikey": ETHANOL_KEY, "smiles": ETHANOL,
         "resolved_from": "PubChem", "lookup": "pubchem:name:SOMETHING ELSE",
         "resolved_title": "ethanol"},
    ])
    rows = resolve(module, PROBES)
    assert [r for r in rows
            if cell(r, "name") == PROBE_WITH_NO_INCHIKEY["name"]] == []


def test_resolve_rejects_a_cached_row_whose_key_disagrees_with_its_smiles(
        tmp_path, no_network):
    # recomputing the InChIKey from the returned SMILES with RDKit and rejecting
    # any that disagrees is mandatory, not optional
    module = module_with_cache(tmp_path / "corrupt", [
        {"name": PROBE_WITH_NO_INCHIKEY["name"], "inchikey": WRONG_KEY,
         "smiles": ETHANOL, "resolved_from": "PubChem",
         "lookup": "pubchem:name:CBK-1234", "resolved_title": "ethanol"},
    ])
    rows = resolve(module, PROBES)
    assert WRONG_KEY not in {cell(r, "inchikey") for r in rows}
    assert rows == []


def test_resolve_records_the_lookup_on_an_accepted_row(tmp_path, no_network):
    module = module_with_cache(tmp_path / "good", [
        {"name": PROBE_WITH_NO_INCHIKEY["name"], "inchikey": ETHANOL_KEY,
         "smiles": ETHANOL, "resolved_from": "PubChem",
         "lookup": "pubchem:name:CBK-1234", "resolved_title": "ethanol"},
    ])
    rows = resolve(module, PROBES)
    assert len(rows) == 1
    row = rows[0]
    assert cell(row, "inchikey") == ETHANOL_KEY
    assert cell(row, "smiles") == ETHANOL
    # a resolved structure did not come from the portal and must never read as
    # if it did, so the lookup is on the row
    values = [text(v) for v in row.values()]
    assert any(re.search(r"pubchem|chembl", v, re.I) for v in values), (
        f"no lookup recorded on the resolved row: {row!r}"
    )


# ----------------------------------------------------------------- write_table

WRITE_COLUMNS = ["inchikey", "target_key", "value", "unit",
                 "assay_description"]
WRITE_ROW = {
    "inchikey": PROBE_WITH_NO_POTENCY_KEY["InChIkey"],
    "target_key": "Q8TEK3",
    "value": "0.06",
    "unit": "nM",
    "assay_description": HAZARD_DESC,
}


@pytest.mark.parametrize("shape", ["rows", "frame"])
def test_write_table_round_trips_a_tab_a_newline_and_a_quote(tmp_path, shape):
    # csv.writer, never '\t'.join: assayDesc carries 5 tabs, 271 line breaks and
    # 10 double quotes, which a naive join spreads across 130 rows
    rows = [dict(WRITE_ROW)]
    data = pd.DataFrame(rows, columns=WRITE_COLUMNS) if shape == "frame" else rows
    preprocess.write_table(tmp_path / f"{preprocess.PREFIX}bioactivity.tsv", data, WRITE_COLUMNS)

    back = read_staged(tmp_path, "bioactivity")
    assert len(back) == 1
    assert back[0]["assay_description"] == HAZARD_DESC
    assert back[0]["value"] == "0.06"


def test_write_table_header_is_the_columns_it_was_given(tmp_path):
    rows = [dict(WRITE_ROW, not_a_column="leaked")]
    preprocess.write_table(tmp_path / f"{preprocess.PREFIX}bioactivity.tsv", rows, WRITE_COLUMNS)
    assert header_of(tmp_path / f"{preprocess.PREFIX}bioactivity.tsv") == WRITE_COLUMNS
    assert "leaked" not in (tmp_path / f"{preprocess.PREFIX}bioactivity.tsv").read_text()


@pytest.mark.parametrize("missing", [None, float("nan")])
def test_write_table_writes_an_empty_cell_for_nothing(tmp_path, missing):
    # a pandas NaN is truthy, so `value or ""` writes the string 'nan'
    rows = [{"inchikey": "X", "unit": missing}]
    preprocess.write_table(tmp_path / f"{preprocess.PREFIX}compound.tsv", rows,
                           ["inchikey", "unit"])
    written = (tmp_path / f"{preprocess.PREFIX}compound.tsv").read_text()
    assert "nan" not in written.lower()
    assert "none" not in written.lower()
    with open(tmp_path / f"{preprocess.PREFIX}compound.tsv", newline="") as handle:
        body = list(csv.reader(handle, delimiter="\t"))[1]
    assert body == ["X", ""]


def test_write_table_writes_an_integer_id_without_a_decimal_point(tmp_path):
    # a column with one missing value is float64 in pandas, which turns
    # canSAR_ID 1354531 into '1354531.0'
    frame = pd.DataFrame({"cansar_id": [1354531, None]})
    preprocess.write_table(tmp_path / f"{preprocess.PREFIX}probe_assessment.tsv", frame, ["cansar_id"])
    with open(tmp_path / f"{preprocess.PREFIX}probe_assessment.tsv", newline="") as handle:
        body = list(csv.reader(handle, delimiter="\t"))[1:]
    assert [row[0] for row in body] == ["1354531", ""]


def test_write_table_writes_the_same_bytes_twice(tmp_path):
    rows = [dict(WRITE_ROW), dict(WRITE_ROW, value="1.5")]
    first, second = tmp_path / "a.tsv", tmp_path / "b.tsv"
    preprocess.write_table(first, rows, WRITE_COLUMNS)
    preprocess.write_table(second, rows, WRITE_COLUMNS)
    assert first.read_bytes() == second.read_bytes()


# --------------------------------------------------------------- write_staging

def test_write_staging_writes_the_contract_and_the_report(tmp_path):
    tables, leftovers = built()
    quarantined = quarantined_numbers()
    preprocess.write_staging(tmp_path, tables, leftovers, quarantined)

    assert header_of(tmp_path / f"{preprocess.PREFIX}compound.tsv") == COMPOUND_COLUMNS
    assert header_of(tmp_path / f"{preprocess.PREFIX}target.tsv") == TARGET_COLUMNS
    assert header_of(tmp_path / f"{preprocess.PREFIX}uniprot.tsv") == UNIPROT_COLUMNS
    assert set(header_of(tmp_path / f"{preprocess.PREFIX}bioactivity.tsv")) == set(BIOACTIVITY_COLUMNS)
    assert header_of(tmp_path / f"{preprocess.PREFIX}probe_assessment.tsv") == ASSESSMENT_COLUMNS
    rejected = read_staged(tmp_path, "rejected_record")
    assert sum(1 for r in rejected if r["stage"] == "bioactivity") == len(quarantined)
    assert isinstance(json.loads((tmp_path / "report.json").read_text()), dict)


def test_main_writes_every_file_of_the_contract(fixture_staging):
    out = fixture_staging()
    for name in preprocess.STAGING_FILES + preprocess.EXTRA_FILES:
        assert (out / f"{preprocess.PREFIX}{name}.tsv").exists(), f"{name}.tsv was not written"
    assert (out / "report.json").exists()


def test_main_validate_prints_the_problems(workdir, fixture_export):
    # the module reports through loguru, like examples/populate_db.py, and loguru
    # writes to the stderr it was imported with, so neither capsys nor capfd sees
    # it. a sink of our own does.
    from loguru import logger

    said = []
    sink = logger.add(said.append, format="{message}")
    out = workdir / "validated"
    out.mkdir(exist_ok=True)
    try:
        with pytest.MonkeyPatch.context() as patch:
            block_network(patch)
            preprocess.main(["--json", str(fixture_export), "--out", str(out),
                             "--validate"])
    finally:
        logger.remove(sink)
    printed = " ".join(said)
    assert "validate:" in printed
    assert "hard errors" in printed
    assert "bioactivity.tsv" in printed          # the counts block a student sees
    assert (out / f"{preprocess.PREFIX}bioactivity.tsv").exists()


def test_pipeline_writes_the_same_bytes_twice(workdir, fixture_export):
    # stable row order, no timestamps in a TSV. report.json is left out of the
    # comparison because a run stamp there would be legitimate.
    directories = []
    for name in ("twice_a", "twice_b"):
        out = workdir / name
        out.mkdir(exist_ok=True)
        with pytest.MonkeyPatch.context() as patch:
            block_network(patch)
            preprocess.main(["--json", str(fixture_export), "--out", str(out)])
        directories.append(out)
    for name in preprocess.STAGING_FILES + preprocess.EXTRA_FILES:
        first = (directories[0] / f"{preprocess.PREFIX}{name}.tsv").read_bytes()
        second = (directories[1] / f"{preprocess.PREFIX}{name}.tsv").read_bytes()
        assert first == second, f"{name}.tsv is not written deterministically"


# ------------------------------------------------------------- staging format

@pytest.mark.parametrize("name,columns", [
    ("compound", COMPOUND_COLUMNS),
    ("target", TARGET_COLUMNS),
    ("uniprot", UNIPROT_COLUMNS),
])
def test_staging_header_is_exactly_the_readme(fixture_staging, name, columns):
    assert header_of(fixture_staging() / f"{preprocess.PREFIX}{name}.tsv") == columns


def test_bioactivity_header_is_the_readme_columns(fixture_staging):
    head = header_of(fixture_staging() / f"{preprocess.PREFIX}bioactivity.tsv")
    assert set(BIOACTIVITY_REQUIRED) <= set(head)
    assert set(head) == set(BIOACTIVITY_COLUMNS)
    assert len(head) == len(set(head))


def test_staging_has_no_column_outside_the_readme(fixture_staging):
    # loader/validate.py reports a column it does not know as 'ignored'
    problems = validate(fixture_staging())
    assert [p for p in problems if "unknown column" in p] == []


def test_target_types_are_in_the_vocabulary(fixture_staging):
    types = {r["type"] for r in read_staged(fixture_staging(), "target")}
    assert types
    assert types <= vocabulary("type")


def test_relations_are_in_the_vocabulary(fixture_staging):
    rows = read_staged(fixture_staging(), "bioactivity")
    written = {r.get("relation") for r in rows} - {"", None}
    assert written <= vocabulary("relation")


def test_chembl_ids_share_one_cell(fixture_staging):
    rows = read_staged(fixture_staging(), "compound")
    cells = {r["inchikey"]: r["chembl_id"] for r in rows}
    assert cells[PROBE_WITH_TWO_CHEMBL_IDS["InChIkey"]].count("|") == 1
    assert all(";" not in value for value in cells.values())


def test_staging_keys_join(fixture_staging):
    out = fixture_staging()
    targets = {r["target_key"] for r in read_staged(out, "target")}
    compounds = {r["inchikey"] for r in read_staged(out, "compound")}
    bioactivity = read_staged(out, "bioactivity")
    assert bioactivity
    assert {r["target_key"] for r in bioactivity} <= targets
    assert {r["inchikey"] for r in bioactivity} <= compounds
    assert {r["target_key"] for r in read_staged(out, "uniprot")} <= targets


def test_a_target_only_under_a_keyless_probe_is_not_written(fixture_staging):
    # 656 accessions in the export, 644 in target.tsv: the difference is the
    # accessions that appear only under a probe with no InChIKey
    targets = {r["target_key"] for r in read_staged(fixture_staging(), "target")}
    assert "P00533" not in targets


# ---------------------------------------------------------------------- report

def test_report_counts_match_the_frames_it_was_given():
    tables, leftovers = built()
    quarantined = quarantined_numbers()
    report = preprocess.report(tables, leftovers, quarantined)
    assert isinstance(report, Mapping)
    for name, frame in tables.items():
        assert count_in(report, name) is not None, (
            f"report() gives no row count for {name}"
        )
        assert count_in(report, name) == len(rows_of(frame))
    assert count_in(report, "rows held for curation") == len(quarantined)


def test_report_accounts_for_every_input_record():
    tables, leftovers = built()
    quarantined = quarantined_numbers()
    report = preprocess.report(tables, leftovers, quarantined)

    compounds = count_in(report, "compound")
    skipped = count_in(report, "compounds skipped")
    written = count_in(report, "bioactivity")
    held = count_in(report, "rows held for curation")
    for name, value in (("compound", compounds), ("skipped_compound", skipped),
                        ("bioactivity", written), ("quarantine", held)):
        assert value is not None, f"report() gives no count for {name}"

    # every probe is either a compound row or a skipped one, once
    assert compounds + skipped == len(PROBES)
    assert compounds == len(KEYED_PROBES)
    assert skipped == 1

    # every number found under a writable compound is either written or held
    # back, and never both (D15)
    assert written + held == len(candidate_numbers())
    assert held == len(quarantined) == 1     # the M-1 s-1 rate constant
    assert written == len(candidate_numbers()) - 1


# ------------------------------------------------------------- end to end

@slow
def test_integration_loads_into_a_fresh_database(portal_db):
    db, report = portal_db()          # loader.load(strict=True) must not raise
    assert report["compounds"] == EXPECTED_COMPOUNDS
    assert db.one("SELECT COUNT(*) FROM compound") == EXPECTED_COMPOUNDS
    assert db.one("SELECT COUNT(*) FROM target") == EXPECTED_TARGETS
    assert db.one("SELECT COUNT(*) FROM uniprot") == EXPECTED_UNIPROTS


@slow
def test_integration_validate_reports_no_hard_error(portal_staging):
    problems = validate_staged(portal_staging())
    hard = [p for p in problems if not p.endswith(("ignored", "kept anyway"))]
    assert hard == []


@slow
def test_integration_foreign_key_check_is_empty(portal_db):
    db, _ = portal_db()
    assert db.conn.execute("PRAGMA foreign_key_check").fetchall() == []


@slow
def test_integration_bioactivity_numbers_are_the_projection(portal_staging,
                                                            portal_db):
    out = portal_staging()
    rows = read_staged(out, "bioactivity")
    # quarantine and skipped_compound are one table now, one row per stage
    rejected = read_staged(out, "rejected_record")
    quarantine = [r for r in rejected if r["stage"] == "bioactivity"]
    # README.md's table calls 3649 the size of bioactivity.tsv while D15
    # holds the 23 quarantined rows back, so the sum is what is asserted here.
    assert len(quarantine) == EXPECTED_QUARANTINED
    assert len(rows) + len(quarantine) == EXPECTED_BIOACTIVITY_NUMBERS

    db, report = portal_db()
    loaded = db.one("SELECT COUNT(*) FROM bioactivity")
    # the loaded count is the file minus what loader/load.py drops, not a number
    # of its own
    assert report["bioactivities"] == len(rows) - EXPECTED_DUPLICATES_DROPPED
    assert loaded == len(rows) - EXPECTED_DUPLICATES_DROPPED


@slow
def test_integration_every_probe_is_written_or_skipped(portal_staging):
    out = portal_staging()
    compounds = read_staged(out, "compound")
    skipped = [r for r in read_staged(out, "rejected_record")
               if r["stage"] == "compound"]
    assert len(compounds) == EXPECTED_COMPOUNDS
    assert len(skipped) == EXPECTED_SKIPPED_COMPOUNDS
    assert len(compounds) + len(skipped) == EXPECTED_PROBES
    keys = [r["inchikey"] for r in compounds]
    assert len(keys) == len(set(keys))


@slow
def test_integration_staging_keys_join(portal_staging):
    out = portal_staging()
    targets = {r["target_key"] for r in read_staged(out, "target")}
    compounds = {r["inchikey"] for r in read_staged(out, "compound")}
    assert {r["target_key"] for r in read_staged(out, "bioactivity")} <= targets
    assert {r["inchikey"] for r in read_staged(out, "bioactivity")} <= compounds
    assert {r["target_key"] for r in read_staged(out, "uniprot")} <= targets


@slow
def test_integration_several_chembl_ids_share_one_cell(portal_staging):
    rows = read_staged(portal_staging(), "compound")
    cells = [r["chembl_id"] for r in rows if r["chembl_id"]]
    assert any("|" in value for value in cells)
    assert all(";" not in value for value in cells)
    ids = [i for value in cells for i in value.split("|")]
    assert all(re.match(r"^CHEMBL[0-9]+$", i) for i in ids)


def test_integration_pinned_reload_reinserts_the_rows_without_a_unit(
        portal_staging):
    out = portal_staging()
    db = ProbeDB(":memory:", create=True)
    load_staging(db, out, strict=True)
    # db.py:297 stores a blank unit as '', load.py:82 looks for None, so the
    # duplicate key of a unitless row never matches what comes back
    stored_without_unit = db.one(
        "SELECT COUNT(*) FROM bioactivity WHERE unit = ''"
    )
    second = load_staging(db, out, strict=True)
    assert second["bioactivities"] == stored_without_unit
    assert second["bioactivities"] < second["duplicates_skipped"]
    # README.md projects 156 rows re-inserted. Measured on the files that
    # are actually written it is 155 unitless rows in bioactivity.tsv, of which
    # 149 reach the database as distinct rows and so come back on every reload.
    rows = read_staged(out, "bioactivity")
    assert len([r for r in rows if not r["unit"]]) == 155
    assert second["bioactivities"] == EXPECTED_REINSERTED_ON_RELOAD
    db.close()


@slow
def test_integration_pinned_duplicate_identity_drops_rows(portal_staging,
                                                             portal_db):
    rows = read_staged(portal_staging(), "bioactivity")
    dropped = len(rows) - len({loader_identity(r) for r in rows})
    # what the drop would be if assay_description and cell_line were in the
    # identity, which is the count of genuine repeats
    genuine = len(rows) - len({
        loader_identity(r) + (r.get("assay_description"), r.get("cell_line"))
        for r in rows
    })
    _, report = portal_db()
    assert report["duplicates_skipped"] == dropped
    assert dropped == EXPECTED_DUPLICATES_DROPPED
    assert genuine == EXPECTED_GENUINE_REPEATS
