"""Chemical Probes Portal export -> a staging directory the loader can read.

    python chemicalprobes.org/preprocess.py \
        --json chemicalprobes.org/ChemicalProbesPortal-6_8_2026.json \
        --out  staging/chemicalprobes.org \
        --db   probe.db

Writes the four files of the staging contract (README section 3), the rows that
need hand curation, the records the schema has no column for, and a report that
reconciles every input record against one of those outputs. With --db it then
loads the staging directory through loader/ into a fresh database.

The mapping this follows, the counts behind it and the decisions are in
preprocess.ipynb and PREPROCESSING.md.

Decided:
  * database/schema.sql is NOT changed, so D2, D3, D5 and D11 are not taken. The
    11 fields with no column, the 1265 in vivo records and the range/error
    numbers do not reach the database. They are still written, beside the staging
    files, so the count is on disk rather than only in a log.
  * database/ and loader/ are NOT touched, so D13 is not taken. Two consequences
    to expect: a reload re-inserts the rows that have no unit, and 35 rows are
    dropped as duplicates of which only 6 are genuine repeats. Both are
    documented in REVIEW_FINDINGS.md and both are asserted by the tests, so they
    are pinned behaviour rather than surprises.
  * the 34 keyless probes are resolved externally (D1). 10 resolve, from the
    cache in resolved_structures.tsv, each re-verified against its SMILES before
    it is used. --resolve refreshes the cache from PubChem.
  * the 260 unsuitable probes get their own entry, unsuitable.tsv, keyed on the
    InChIKey.
"""

import argparse
import csv
import json
import re
import sys
from pathlib import Path

import pandas as pd
from loguru import logger

SOURCE_DB = "Chemical Probes Portal"
PORTAL_PREFIX = "https://www.chemicalprobes.org/"
DOI_PREFIX = "https://doi.org/"
PUBMED_PREFIX = "https://pubmed.ncbi.nlm.nih.gov/"

# every portal target entry is one accession, so a one-member protein target
TARGET_TYPE = "protein"

STAGING_FILES = ("compound", "target", "uniprot", "bioactivity")
# unsuitable is the one that carries a foreign key: its inchikey references
# compound, the way chembl does
EXTRA_FILES = ("unsuitable", "quarantine", "compound_annotation",
               "target_annotation", "in_vivo", "reference", "skipped_compound")

RESOLVED_CACHE = Path(__file__).resolve().parent / "resolved_structures.tsv"

# the 18 probe keys the mapping in preprocess.ipynb covers. a new one has to
# stop the run, not be ignored
PROBE_KEYS = (
    "name", "rating_in_cell", "rating_in_organism", "rating_count", "unsuitable",
    "URL", "primary_targets", "published_date", "InChIkey", "smiles", "pains",
    "toxicophore", "canSAR_ID", "ChEMBL_ID", "PMID", "inVitroValidations",
    "inVivoValidations", "control_compounds",
)
TIERS = (("in vitro", "inVitroValidations"), ("in cell", "inCellValidations"))

COLUMNS = {
    "compound": ["inchikey", "smiles", "chembl_id", "name"],
    "target": ["target_key", "type", "name"],
    "uniprot": ["uniprot_id", "target_key", "hgnc", "species"],
    "bioactivity": ["inchikey", "target_key", "moa", "bioactivity_type", "relation",
                    "value", "unit", "assay_type", "assay_description", "cell_line",
                    "concentration", "concentration_unit", "source_db", "source",
                    "source_xref", "xref_id"],
    "unsuitable": ["inchikey", "name", "smiles", "chembl_id", "cansar_id",
                   "portal_path", "published_date", "pains", "toxicophore",
                   "rating_in_cell", "rating_in_organism", "rating_count",
                   "reference", "source_db"],
    "quarantine": ["inchikey", "target_key", "reason", "raw", "fragment", "relation",
                   "value", "unit", "bioactivity_type", "assay_description",
                   "source_db", "source"],
    "compound_annotation": ["inchikey", "source_db", "property", "ordinal", "value"],
    "target_annotation": ["target_key", "source_db", "property", "ordinal", "value"],
    "in_vivo": ["inchikey", "organism", "dose_value", "dose_unit", "route", "dose_raw",
                "source_db", "source"],
    "reference": ["inchikey", "xref_id", "source_xref", "raw"],
    "skipped_compound": ["name", "reason", "portal_path", "targets", "validations"],
}


# ---------------------------------------------------------------- reading

def check_vocabularies():
    """The closed vocabularies live in the CHECK constraints in schema.sql and
    nowhere else, so what this module writes has to be read back out of them."""
    from probedb.schema import vocabulary

    for name, used in (("relation", set(RELATION.values())),
                       ("type", {TARGET_TYPE})):
        allowed = vocabulary(name)
        if not used <= allowed:
            raise ValueError(
                f"{name}: this module writes {sorted(used - allowed)}, which "
                f"schema.sql does not allow ({sorted(allowed)})"
            )


def load_export(path):
    payload = json.loads(Path(path).read_text())
    if not isinstance(payload, dict) or "probes" not in payload:
        raise ValueError(f"{path}: no 'probes' key, is this a portal export?")
    probes = payload["probes"]
    if not isinstance(probes, list):
        raise ValueError(f"{path}: 'probes' is {type(probes).__name__}, not a list")
    for i, probe in enumerate(probes):
        if not isinstance(probe, dict):
            raise ValueError(f"{path}: probe {i} is {type(probe).__name__}, not an object")
        missing = set(PROBE_KEYS) - set(probe)
        if missing:
            raise ValueError(f"{path}: probe {i} is missing {sorted(missing)}")
    unknown = {key for probe in probes for key in probe} - set(PROBE_KEYS)
    if unknown:
        raise ValueError(
            f"{path}: unmapped probe key(s) {sorted(unknown)}. A new key means the "
            "mapping in preprocess.ipynb is out of date"
        )
    return probes


def flatten(probes):
    frame = {}
    frame["probe"] = pd.DataFrame([
        {"probe_ix": i,
         "name": as_text(p["name"]),
         "inchikey": as_text(p["InChIkey"]).upper(),
         "smiles": as_text(p["smiles"]),
         "url": as_text(p["URL"]),
         "published_date": as_text(p["published_date"]),
         "unsuitable": as_text(p["unsuitable"]),
         "pains": as_text(p["pains"]),
         "toxicophore": as_text(p["toxicophore"]),
         "cansar_id": as_number(p["canSAR_ID"]),
         "rating_in_cell": as_number(p["rating_in_cell"]),
         "rating_in_organism": as_number(p["rating_in_organism"]),
         "rating_count": as_number(p["rating_count"]),
         "structure_from": "",
         "structure_lookup": ""}
        for i, p in enumerate(probes)
    ], columns=["probe_ix", "name", "inchikey", "smiles", "url", "published_date",
                "unsuitable", "pains", "toxicophore", "cansar_id", "rating_in_cell",
                "rating_in_organism", "rating_count", "structure_from",
                "structure_lookup"])
    frame["target"] = pd.DataFrame([
        {"probe_ix": i, "target_ix": j,
         "symbol": as_text(t["name"]),
         "uniprot_id": as_text(t["uniprot_id"]).upper(),
         "target_class": as_text(t["class"]),
         "subclass": as_text(t["subClass"]),
         "moa": normalise_moa(t["moa"])}
        for i, p in enumerate(probes) for j, t in enumerate((p.get("primary_targets") or []))
    ], columns=["probe_ix", "target_ix", "symbol", "uniprot_id", "target_class",
                "subclass", "moa"])
    frame["validation"] = pd.DataFrame([
        {"probe_ix": i, "target_ix": j, "tier": tier,
         "potency": v.get("potency"),
         "potency_value": v.get("potencyValue"),
         "assay_desc": v.get("assayDesc")}
        for i, p in enumerate(probes) for j, t in enumerate((p.get("primary_targets") or []))
        for tier, key in TIERS for v in (t.get(key) or [])
    ], columns=["probe_ix", "target_ix", "tier", "potency", "potency_value",
                "assay_desc"])
    frame["invivo"] = pd.DataFrame([
        # 7 records carry no 'dose' key at all, so .get and not ["dose"]
        {"probe_ix": i, "organism": as_text(v.get("organism")),
         "dose": as_text(v.get("dose"))}
        for i, p in enumerate(probes) for v in (p.get("inVivoValidations") or [])
    ], columns=["probe_ix", "organism", "dose"])
    frame["chembl"] = pd.DataFrame([
        {"probe_ix": i, "chembl_id": as_text(c).upper()}
        for i, p in enumerate(probes) for c in (p.get("ChEMBL_ID") or []) if as_text(c)
    ], columns=["probe_ix", "chembl_id"])
    frame["reference"] = pd.DataFrame([
        {"probe_ix": i, "ref": as_text(r)}
        for i, p in enumerate(probes) for r in (p.get("PMID") or []) if as_text(r)
    ], columns=["probe_ix", "ref"])
    frame["control"] = pd.DataFrame([
        {"probe_ix": i, "control_name": as_text(c)}
        for i, p in enumerate(probes) for c in (p.get("control_compounds") or []) if as_text(c)
    ], columns=["probe_ix", "control_name"])
    return frame


# ------------------------------------------------------- normalising values

def as_text(value):
    return "" if value is None or value != value else str(value).strip()


def as_number(value):
    """A number as it was written: 0 not 0.0, 1354531 not 1354531.0."""
    if value is None or value != value or value == "":
        return ""
    try:
        number = float(value)
    except (TypeError, ValueError):
        return as_text(value)
    return str(int(number)) if number.is_integer() else str(number)


MOA_FIX = {
    "antagoinist, degrader": "antagonist, degrader",   # typo in the export
    "molecular glues": "molecular glue",               # plural of one mechanism
}


def normalise_moa(value):
    moa = re.sub(r"\s+", " ", as_text(value)).lower()
    return MOA_FIX.get(moa, moa)


ENDPOINT_FIX = {
    "ic 50": "IC50", "ec 50": "EC50", "dc 50": "DC50", "gi 50": "GI50",
    "gi50*": "GI50", "k 0.5": "K0.5", "kd apparent": "Kd(app)", "dc50 @ 16h": "DC50",
    "inh": "% inhibition", "inhibition": "% inhibition", "% in": "% inhibition",
    "% ac": "% activation", "pkbb": "pKb", "not done": None, "delta tm": "ΔTm",
    "activity": "activity", "residual activity": "residual activity",
    "ratio": "ratio",
}
CANON_ENDPOINT = {e.lower(): e for e in
                  ["IC50", "EC50", "DC50", "GI50", "Ki", "Kd", "Dmax", "ΔTm", "pIC50",
                   "pEC50", "pGI50", "pKb", "IC30", "K0.5", "Kd(app)", "MEC1.5",
                   "MEC2.0", "Emax"]}


def normalise_endpoint(value):
    endpoint = re.sub(r"\s+", " ", as_text(value))
    if not endpoint or endpoint.lower() in ("none", "nan"):
        return None
    low = endpoint.lower()
    if low in ENDPOINT_FIX:
        return ENDPOINT_FIX[low]
    if low in CANON_ENDPOINT:
        return CANON_ENDPOINT[low]
    return re.sub(r"\b(IC|EC|DC|GI) ?50\b", r"\g<1>50", endpoint, flags=re.I)


def split_endpoints(value):
    return [e for e in (normalise_endpoint(p) for p in re.split(r"[,;]", as_text(value)))
            if e]


def assign_endpoint(fragment, potency):
    labels = split_endpoints(potency)
    if len(labels) <= 1:
        return labels[0] if labels else None
    for label in labels:
        loose = re.escape(label).replace(r"\ ", r"\s*").replace("50", r"\s*50")
        if re.search(loose, fragment, re.I):
            return label
    return None


# rate constants and reciprocal units first, they contain the shorter units. the
# ambiguous single letters M, K and h must not sit inside a word or follow '(',
# or the M of 'T790M' is molar, a '(K)' cell line is a temperature and an '(H)'
# species label is an hour
UNIT_ALT = [
    (r"M\s*[−–-]\s*1\s*[·.]?\s*[sS]\s*[−–-]\s*1"
     r"|[sS]\s*[−–-]\s*1\s*[·.]?\s*M\s*[−–-]\s*1"
     r"|/\(mol/L\)\s*s|per\s*s/[uµμ]mol/L|per\s*M\s*per\s*sec", "M-1s-1"),
    (r"min\s*[−–-]\s*1", "min-1"),
    (r"nmol/L", "nM"), (r"[uµμ]mol/L", "uM"), (r"mmol/L", "mM"), (r"mol/L", "M"),
    (r"pM\b", "pM"), (r"nM\b", "nM"), (r"[uµμ]M\b", "uM"), (r"mM\b", "mM"),
    (r"(?<![A-Za-z0-9(])M\b", "M"),
    (r"%|Percent\s*of\s*Control", "%"),
    (r"°\s*C|ºC|Celciuys|Celsius|degrees|Kelvin|(?<![A-Za-z0-9(])K\b", "degC"),
    (r"fold|(?<=\d)x\b", "fold"),
    (r"(?<![A-Za-z0-9(])h\b|hours?\b", "h"),
    (r"(?<![A-Za-z0-9(])min\b", "min"),
]
UNIT_RE = re.compile("|".join(f"(?P<u{i}>{p})" for i, (p, _) in enumerate(UNIT_ALT)), re.I)
UNIT_CANON = {f"u{i}": c for i, (_, c) in enumerate(UNIT_ALT)}
ANY_UNIT = "|".join(p for p, _ in UNIT_ALT)
MOLAR = {"pM", "nM", "uM", "mM", "M"}

NUM = r"\d+(?:\.\d+)?"
RELATION = {"<": "<", ">": ">", "~": "~", "=": "=", "≤": "<=", "≥": ">=",
            "<=": "<=", ">=": ">=", "≈": "~"}
WORD_RELATION = [                       # a censored value must not become exact
    (r"\b(?:less than|below|under|up to|at most|fewer than)\b", "<"),
    (r"\b(?:greater than|more than|above|over|at least)\b", ">"),
    (r"\b(?:about|around|approximately|approx\.?|circa|ca\.)\b", "~"),
]
NOT_A_VALUE = {"", "none", "na", "n/a", "nd", "not determined", "not done",
               "-", "not available"}
P_SCALE = re.compile(r"^p[A-Z]", re.I)
PERCENT_ENDPOINT = re.compile(r"^(?:%|dmax|emax)", re.I)
# only a real concentration is lifted out: '88% @ 20 h' is a time window, and
# 20 hours in bioactivity.concentration would be a lie
CONCENTRATION_UNIT = "|".join(
    p for p, canon in UNIT_ALT if canon in ("pM", "nM", "uM", "mM", "M", "%"))
AT_CONCENTRATION = re.compile(
    rf"(?:@|\bat\b)\s*(?P<value>{NUM})\s*(?P<unit>{CONCENTRATION_UNIT})", re.I)
EXPONENT = re.compile(r"[xX×*]\s*10\s*[\^*]?\s*[-−–+]?\d+|[eE][-+]\d+")
SLASHED_NUMBERS = re.compile(rf"(?<![A-Za-z0-9]){NUM}\s*/\s*{NUM}(?![A-Za-z0-9])")
FRAGMENT = re.compile(
    rf"""(?P<rel>[<>~≤≥≈]=?|=)?\s*
         (?<![A-Za-z0-9.])(?P<lo>{NUM})
         (?:\s*(?:±|\+/-|\+-)\s*(?P<err>{NUM}))?
         (?:\s*(?:{ANY_UNIT})?\s*(?:[-–—]|\bto\b)\s*(?P<hi>{NUM})(?![0-9.]))?""",
    re.X | re.I,
)


def canonical_unit(text):
    match = UNIT_RE.search(as_text(text))
    return UNIT_CANON[match.lastgroup] if match else None


def repair(text):
    text = re.sub(r"(?<=\d)\.\s+(?=\d)", ".", text)          # '0. 69'  -> '0.69'
    text = re.sub(r"(?<=\d),(?=\d{3}\b)", "", text)          # '1,100'  -> '1100'
    text = re.sub(r"(?<=\d),(?=\d(?!\d{2}))", ".", text)     # '6,5 nM' -> '6.5 nM'
    text = re.sub(r"(?<=\d)\s*n\s+M\b", " nM", text)         # '85n M'  -> '85 nM'
    for pattern, symbol in WORD_RELATION:
        text = re.sub(pattern, symbol, text, flags=re.I)
    return text


def split_value_fragments(text):
    out = []
    for part in (p for p in re.split(r";(?!\d)", text) if p.strip()):   # keep 'MV4;11'
        pieces = [x for x in re.split(r",", part) if x.strip()]
        if len(pieces) == 1:
            # '760nM and 1000nM', '1 or 41 nM', '7.8 nM (human) 2.4 nM (Mouse)'
            pieces = [x for x in re.split(r"\b(?:and|or)\b|(?<=[)\]])\s+(?=[<>~]?\s*\d)",
                                          part, flags=re.I) if x and x.strip()]
            if len(pieces) > 1 and not all(re.search(r"\d", x) for x in pieces):
                pieces = [part]
        out.extend(pieces)
    return out


def quarantine_reason(row, potency, exponent, slashed):
    if exponent:
        return "scientific-notation factor not applied"
    if row["unit"] in ("M-1s-1", "min-1"):
        return "rate constant, not a potency"
    if slashed:
        return "two numbers separated by '/'"
    if row["unit"] in ("h", "min"):
        return "a duration, not a potency"
    if row["value"] == 0 and row["unit"] in MOLAR:
        return "zero potency"
    if row["unit"] == "M" and row["value"] > 1:
        return "implausible molar value"
    if row["value_high"] is not None and row["value_high"] < row["value"]:
        return "range reads high to low"
    if P_SCALE.match(as_text(potency)) and row["unit"] in MOLAR:
        return "p-scale endpoint with a molar unit"
    return None


def parse_potency_value(raw, potency=None):
    text = as_text(raw)
    if text.lower() in NOT_A_VALUE:
        return []
    exponent, slashed = bool(EXPONENT.search(text)), bool(SLASHED_NUMBERS.search(text))
    text = repair(text)
    rows = []
    for fragment in split_value_fragments(text):
        # the concentration an assay ran at is not the measurement
        at = AT_CONCENTRATION.search(fragment)
        measured = fragment[:at.start()] + fragment[at.end():] if at else fragment
        match = FRAGMENT.search(measured)
        if not match:
            continue
        rows.append({
            "relation": RELATION.get(match.group("rel") or "", None),
            "value": float(match.group("lo")),
            "value_high": float(match.group("hi")) if match.group("hi") else None,
            "error": float(match.group("err")) if match.group("err") else None,
            "unit": canonical_unit(measured[match.end():]) or canonical_unit(measured),
            "unit_inherited": False,
            "concentration": float(at.group("value")) if at else None,
            "concentration_unit": canonical_unit(at.group("unit")) if at else None,
            "bioactivity_type": assign_endpoint(fragment, potency),
            "fragment": fragment.strip(),
            "raw": as_text(raw),
            "quarantine": None,
        })
    if not rows:
        return []

    labels = split_endpoints(potency)
    single = len(labels) <= 1
    trailing = [r["unit"] for r in rows if r["unit"]]
    if trailing and single:                            # 'BD1 85.5, BD2 220 nM'
        for row in rows:
            if not row["unit"]:
                row["unit"], row["unit_inherited"] = trailing[-1], True
    if not trailing and single:
        label = labels[0] if labels else ""
        fallback = ("-log(M)" if P_SCALE.match(label)
                    else "%" if PERCENT_ENDPOINT.match(label) else None)
        if fallback:
            for row in rows:
                row["unit"], row["unit_inherited"] = fallback, True
    for row in rows:
        row["quarantine"] = quarantine_reason(row, potency, exponent, slashed)
    return rows


DOSE = re.compile(
    r"(?P<value>\d+(?:\.\d+)?)\s*(?P<unit>mg/[Kk]g(?:/[Dd]ay)?|mgKg|mpk|mg/m[Ll]|"
    r"[uµμ]M/[Kk]g|[uµμ]mol/[Kk]g|nM/kg|mg\b|[uµμ]g\b)", re.I)
ROUTE = re.compile(r"\b(IV|PO|IP|SC|oral|gavage|topical)\b", re.I)
DOSE_UNIT = {"mgkg": "mg/kg", "mpk": "mg/kg", "mg/kg/day": "mg/kg/day",
             "mg/ml": "mg/mL", "um/kg": "uM/kg", "umol/kg": "umol/kg",
             "nm/kg": "nM/kg", "mg": "mg", "ug": "ug"}


def normalise_dose_unit(unit):
    key = as_text(unit).lower().replace("µ", "u").replace("μ", "u")
    return DOSE_UNIT.get(key, key)


def parse_dose(raw):
    text = as_text(raw)
    if not text or text.lower() in NOT_A_VALUE or text.lower() == "unknown":
        return []
    matches = list(DOSE.finditer(text))
    out = []
    for i, match in enumerate(matches):
        # the route follows its own dose. looking backwards would pick up the
        # route of the dose before it: '1 mg/Kg IV, 5 mg/Kg' has one route
        end = matches[i + 1].start() if i + 1 < len(matches) else len(text)
        route = ROUTE.search(text[match.end():end])
        value = float(match.group("value"))
        # a dose range keeps its low end, the same convention bioactivity uses
        low = re.search(rf"({NUM})\s*(?:[-–—]|\bto\b)\s*{re.escape(match.group('value'))}\b",
                        text[:match.end()])
        if low:
            value = float(low.group(1))
        out.append((value, normalise_dose_unit(match.group("unit")),
                    normalise_route(route.group(1)) if route else None))
    return out


ROUTE_CANON = {"oral": "PO", "gavage": "PO"}


def normalise_route(route):
    route = as_text(route).lower()
    return ROUTE_CANON.get(route, route.upper())


DOI = re.compile(r"(10\.\d{4,9}/\S+)")
PMID_NUM = re.compile(r"pubmed[^\d]*?(\d{6,9})")


def reference_source(value):
    text = as_text(value)
    doi = DOI.search(text)
    if doi:
        return DOI_PREFIX, doi.group(1).rstrip(").;,")
    pmid = PMID_NUM.search(text)
    if pmid:
        return PUBMED_PREFIX, pmid.group(1)
    return None, text


def paper_label(value):
    """A reference -> the label the template puts in bioactivity_source.source."""
    prefix, xref = reference_source(value)
    if prefix == DOI_PREFIX:
        return "DOI:" + xref
    if prefix == PUBMED_PREFIX:
        return "PMID:" + xref
    return None


def portal_source(url):
    text = as_text(url)
    if not text.startswith(PORTAL_PREFIX):
        return None, text
    # the whole path, not the last segment: 260 probes live under /unsuitables/
    return PORTAL_PREFIX, text[len(PORTAL_PREFIX):]


def cited_paper(refs):
    """The one paper a probe cites, or None when it cites several or none.

    A publisher url with no DOI in it counts as a paper of its own, because it
    might be a second one. Only an unambiguous single paper earns the attribution.
    """
    papers = {paper_label(r) or as_text(r).lower() for r in refs if as_text(r)}
    if len(papers) != 1:
        return None
    only = papers.pop()
    return only if only.startswith(("DOI:", "PMID:")) else None

# ---------------------------------------------------------------- resolving

def read_resolved_cache(path=RESOLVED_CACHE):
    """The externally resolved structures, only those that verify (D1).

    A cached row is used only if the InChIKey recomputed from its SMILES matches
    the one cached with it. Without RDKit nothing verifies, so nothing is used --
    an unverified structure is worse than a missing one.
    """
    path = Path(path)
    if not path.exists():
        logger.warning(f"no {path.name}, the probes with no InChIKey stay skipped")
        return {}
    try:
        from rdkit import Chem, RDLogger
        from rdkit.Chem.inchi import MolToInchiKey

        RDLogger.DisableLog("rdApp.*")
    except ImportError:
        logger.warning(f"rdkit not installed, ignoring {path.name}")
        return {}

    resolved = {}
    with open(path, newline="") as handle:
        for row in csv.DictReader(handle, delimiter="\t"):
            name, key, smiles = (as_text(row.get(c)) for c in ("name", "inchikey", "smiles"))
            if not (name and key and smiles):
                continue
            mol = Chem.MolFromSmiles(smiles)
            if mol is None or MolToInchiKey(mol) != key:
                logger.warning(
                    f"{path.name}: {name} smiles and inchikey disagree, ignored")
                continue
            resolved[name] = {
                "inchikey": key.upper(),
                "smiles": smiles,
                "resolved_from": as_text(row.get("resolved_from")) or "unknown",
                "lookup": as_text(row.get("lookup")),
            }
    return resolved


def resolve_missing_structures(probe, cache=RESOLVED_CACHE):
    """The structures the export leaves empty, filled in from the cache (D1).

    Returns one row per probe that resolved -- name, inchikey, smiles and the
    lookup it came from -- and nothing for the ones that did not. Empty when
    there is no cache, so the pipeline never depends on the network. A resolved
    structure must never read as a portal one, which is what resolved_from and
    lookup record.
    """
    resolved = read_resolved_cache(cache)
    if not resolved:
        return []
    # a portal structure is never overwritten, only a missing one filled in
    return [{"name": as_text(row["name"]), **resolved[as_text(row["name"])]}
            for row in probe.to_dict("records")
            if not as_text(row["inchikey"]) and as_text(row["name"]) in resolved]


def apply_resolved(probe, resolved):
    """Put the resolved structures onto the probe frame."""
    by_name = {row["name"]: row for row in resolved}
    probe = probe.copy()
    for row in probe[probe.inchikey == ""].itertuples():
        hit = by_name.get(row.name)
        if not hit:
            continue
        probe.loc[row.Index,
                  ["inchikey", "smiles", "structure_from", "structure_lookup"]] = [
            hit["inchikey"], hit["smiles"], hit["resolved_from"], hit["lookup"]]
    return probe


def refresh_resolved_cache(probes, path=RESOLVED_CACHE):
    """Re-query PubChem, then ChEMBL, for the probes with no structure (--resolve).

    Only writes rows whose SMILES and InChIKey agree, so the cache is verified
    when it is written as well as when it is read.
    """
    import time
    import urllib.parse
    import urllib.request

    from rdkit import Chem, RDLogger
    from rdkit.Chem.inchi import MolToInchiKey

    RDLogger.DisableLog("rdApp.*")

    def fetch(url):
        try:
            with urllib.request.urlopen(url, timeout=20) as response:
                return json.load(response)
        except Exception:
            return None

    rows = []
    for probe in probes:
        if as_text(probe["InChIkey"]):
            continue
        name = as_text(probe["name"])
        payload = fetch("https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/name/"
                        + urllib.parse.quote(name)
                        + "/property/SMILES,InChIKey,Title/JSON")
        properties = (payload or {}).get("PropertyTable", {}).get("Properties", [])
        hit = properties[0] if properties else None
        source, lookup = "PubChem", f"pubchem:name:{name}"
        if not hit and probe["ChEMBL_ID"]:
            chembl_id = as_text(probe["ChEMBL_ID"][0])
            payload = fetch(f"https://www.ebi.ac.uk/chembl/api/data/molecule/"
                            f"{chembl_id}?format=json")
            structures = (payload or {}).get("molecule_structures") or {}
            hit = ({"InChIKey": structures.get("standard_inchi_key"),
                    "SMILES": structures.get("canonical_smiles"),
                    "Title": (payload or {}).get("pref_name")}
                   if structures.get("standard_inchi_key") else None)
            source, lookup = "ChEMBL", f"chembl:id:{chembl_id}"
        time.sleep(0.25)
        if not (hit and hit.get("InChIKey") and hit.get("SMILES")):
            continue
        mol = Chem.MolFromSmiles(hit["SMILES"])
        if mol is None or MolToInchiKey(mol) != hit["InChIKey"]:
            continue
        rows.append({"name": name, "inchikey": hit["InChIKey"], "smiles": hit["SMILES"],
                     "resolved_from": source, "lookup": lookup,
                     "resolved_title": as_text(hit.get("Title"))})
    if not rows:
        # the network was down, or every lookup failed. an empty cache is worse
        # than the one already on disk, so the file is left alone
        raise RuntimeError(f"--resolve resolved nothing, {Path(path).name} untouched")
    write_table(Path(path), rows,
                ["name", "inchikey", "smiles", "resolved_from", "lookup", "resolved_title"])
    return rows


# ------------------------------------------------------- building the tables

def build_compound(probe, chembl):
    ids = (chembl.groupby("probe_ix").chembl_id
           .apply(lambda s: "|".join(dict.fromkeys(s))))   # loader splits on '|'
    rows, seen = [], set()
    for row in probe[probe.inchikey != ""].itertuples():
        if row.inchikey in seen:
            # two probes on one key: the loader would keep the first row's name
            # and silently attribute the second probe's measurements to it
            logger.warning(f"{row.name}: {row.inchikey} is already used, not written")
            continue
        seen.add(row.inchikey)
        rows.append({"inchikey": row.inchikey, "smiles": row.smiles,
                     "chembl_id": ids.get(row.probe_ix, ""), "name": row.name})
    return rows


def build_target(target):
    """One row per distinct accession. type is always 'protein': every portal
    target entry is a single accession, so it is a one-member target.

    The caller passes only the entries whose probe is being written.
    """
    rows = {}
    for row in target.itertuples():
        rows.setdefault(row.uniprot_id,
                        {"target_key": row.uniprot_id, "type": TARGET_TYPE,
                         "name": row.symbol})
    return [rows[key] for key in sorted(rows)]


def build_uniprot(target):
    """species and entrez_gene stay empty: neither is in the export (D12)."""
    rows = {}
    for row in target.itertuples():
        rows.setdefault(row.uniprot_id,
                        {"uniprot_id": row.uniprot_id, "target_key": row.uniprot_id,
                         "hgnc": row.symbol, "species": ""})
    return [rows[key] for key in sorted(rows)]


BINDING = re.compile(r"\b(?:SPR|ITC|BROMOscan|thermal shift|DSF|MST|NanoBRET|"
                     r"fluorescence polarization|radioligand|scintillation)\b", re.I)


def assay_type_for(tier, description):
    """The portal's 'in vitro' tier means cell-free, which is 'biochemical' here.

    Refined to 'binding' where the description names a binding technique, which
    is the distinction the template's assay_type makes (D7).
    """
    if tier == "in cell":
        return "cell"
    return "binding" if BINDING.search(as_text(description)) else "biochemical"


def build_provenance(probe, reference):
    """probe_ix -> (source, source_xref, xref_id), following the template.

    The paper goes in source when the probe cites exactly one, because that is
    where staging/_template puts an aggregator's reference. Otherwise the probe
    names itself, since picking one of several cited papers would invent an
    attribution. source_xref is always the portal record the number was read
    from, so source_url resolves either way.
    """
    refs = ({} if reference is None or len(reference) == 0
            else reference.groupby("probe_ix").ref.apply(list))
    out = {}
    for row in probe.itertuples():
        prefix, path = portal_source(row.url)
        out[row.probe_ix] = (cited_paper(refs.get(row.probe_ix, [])) or row.name,
                             path, prefix or "")
    return out


def build_bioactivity(probe, target, validation, reference=None):
    """One row per parsed number that can be trusted (D15)."""
    return split_bioactivity(probe, target, validation, reference)[0]


def split_bioactivity(probe, target, validation, reference=None):
    """(loadable rows, held-back rows). One row per parsed number."""
    provenance = build_provenance(probe, reference)
    # a Series attribute called .name is the index label in pandas, so the probe
    # frame is read through plain dicts here and nowhere through .loc[...].name
    key_of = {r.probe_ix: r.inchikey for r in probe.itertuples() if r.inchikey}
    accession = {(r.probe_ix, r.target_ix): r.uniprot_id for r in target.itertuples()}
    moa = {(r.probe_ix, r.target_ix): r.moa for r in target.itertuples()}
    rows, quarantined = [], []
    for record in validation.itertuples():
        if record.probe_ix not in key_of:
            continue
        source, source_xref, xref_id = provenance[record.probe_ix]
        fragments = parse_potency_value(record.potency_value, record.potency)
        description = as_text(record.assay_desc)
        for fragment in fragments or [None]:
            row = bioactivity_row_from(
                fragment, description,
                inchikey=key_of[record.probe_ix],
                target_key=accession[(record.probe_ix, record.target_ix)],
                moa=moa[(record.probe_ix, record.target_ix)],
                assay_type=assay_type_for(record.tier, description),
                provenance=(source, source_xref, xref_id),
                siblings=len(fragments),
            )
            if fragment and fragment["quarantine"]:
                quarantined.append({
                    "inchikey": row["inchikey"], "target_key": row["target_key"],
                    "reason": fragment["quarantine"], "raw": fragment["raw"],
                    "fragment": fragment["fragment"], "relation": row["relation"],
                    "value": row["value"], "unit": row["unit"],
                    "bioactivity_type": row["bioactivity_type"],
                    "assay_description": row["assay_description"],
                    "source_db": SOURCE_DB,
                    "source": row["source"]})
                continue
            rows.append(row)
    return rows, quarantined


def bioactivity_row_from(fragment, description="", inchikey="", target_key="",
                         moa="", assay_type="", provenance=("", "", ""),
                         siblings=1):
    """One parsed number -> one bioactivity row.

    The schema has no column for the upper end of a range or for a standard
    deviation, so the fragment text is appended to the description whenever it
    carries either -- otherwise a 10-100 nM range would be written as a bare
    10 nM and read as exact. A range is also censored with '>=', which is a true
    statement about it, where no relation at all is not.
    """
    fragment = fragment or {}
    source, source_xref, xref_id = provenance
    spread = fragment.get("value_high") is not None or fragment.get("error") is not None
    note = description
    label = as_text(fragment.get("fragment"))
    if label and label != description and (siblings > 1 or spread):
        note = f"{description} | {label}".strip(" |")
    relation = fragment.get("relation") or ""
    if not relation and fragment.get("value_high") is not None:
        relation = ">="
    return {
        "inchikey": inchikey,
        "target_key": target_key,
        "moa": moa,
        "bioactivity_type": fragment.get("bioactivity_type") or "",
        "relation": relation,
        "value": as_number(fragment.get("value")),
        "unit": fragment.get("unit") or "",
        "assay_type": assay_type,
        "assay_description": note,
        "cell_line": "",                                  # D10
        "concentration": as_number(fragment.get("concentration")),
        "concentration_unit": fragment.get("concentration_unit") or "",
        "source_db": SOURCE_DB,
        "source": source,
        "source_xref": source_xref,
        "xref_id": xref_id,
    }


def build_unsuitable(probe, chembl, reference):
    """The whole portal record for the compounds it has ruled out.

    260 rows keyed on inchikey, which references compound. They carry no target,
    no validation and no in vivo record, so a flat entry loses nothing.
    """
    ids = chembl.groupby("probe_ix").chembl_id.apply(lambda s: "|".join(dict.fromkeys(s)))
    refs = reference.groupby("probe_ix").ref.apply(
        lambda s: "|".join(reference_url(v) for v in s))
    rows = []
    for row in probe[(probe.unsuitable == "Yes") & (probe.inchikey != "")].itertuples():
        prefix, path = portal_source(row.url)
        rows.append({"inchikey": row.inchikey, "name": row.name, "smiles": row.smiles,
                     "chembl_id": ids.get(row.probe_ix, ""),
                     "cansar_id": row.cansar_id, "portal_path": path,
                     "published_date": row.published_date, "pains": row.pains,
                     "toxicophore": row.toxicophore,
                     "rating_in_cell": row.rating_in_cell,
                     "rating_in_organism": row.rating_in_organism,
                     "rating_count": row.rating_count,
                     "reference": refs.get(row.probe_ix, ""), "source_db": SOURCE_DB})
    return rows


def reference_url(value):
    prefix, xref = reference_source(value)
    return (prefix + xref) if prefix else xref

# ------------------------------------------------------------ what is left

def build_leftovers(probe, target, invivo, reference, control, validation=None):
    """The records the schema has no column for, one list per file.

    Written whether or not the schema ever grows the tables, so the count is on
    disk rather than only in a log (D2, D3, D4, D5).
    """
    keyed = set(probe.loc[probe.inchikey != "", "probe_ix"])
    key_of = dict(zip(probe.probe_ix, probe.inchikey))
    name_of = dict(zip(probe.probe_ix, probe.name))
    path_of = {r.probe_ix: portal_source(r.url)[1] for r in probe.itertuples()}

    annotation = []
    for row in probe[probe.inchikey != ""].itertuples():
        for prop, value in (("rating_in_cell", row.rating_in_cell),
                            ("rating_in_organism", row.rating_in_organism),
                            ("rating_count", row.rating_count),
                            ("unsuitable", row.unsuitable),
                            ("pains", row.pains),
                            ("toxicophore", row.toxicophore),
                            ("published_date", row.published_date),
                            ("canSAR_ID", row.cansar_id),
                            ("url", row.url)):
            annotation.append({"inchikey": row.inchikey, "source_db": SOURCE_DB,
                               "property": prop, "ordinal": 0, "value": value})
        # a resolved structure did not come from the portal and must say so
        if row.structure_from:
            annotation.append({"inchikey": row.inchikey, "source_db": row.structure_from,
                               "property": "structure_source", "ordinal": 0,
                               "value": row.structure_lookup})
    for probe_ix, group in control[control.probe_ix.isin(keyed)].groupby("probe_ix"):
        for ordinal, name in enumerate(group.control_name):
            annotation.append({"inchikey": key_of[probe_ix], "source_db": SOURCE_DB,
                               "property": "control_compound", "ordinal": ordinal,
                               "value": name})

    # the taxonomy belongs to the portal's (probe, target) row, so one accession
    # carries several values: keep every distinct one, numbered, rather than
    # letting the last row seen win
    distinct = {}
    for row in target[target.probe_ix.isin(keyed)].itertuples():
        for prop, value in (("class", row.target_class), ("subClass", row.subclass)):
            if value:
                distinct.setdefault((row.uniprot_id, prop), [])
                if value not in distinct[(row.uniprot_id, prop)]:
                    distinct[(row.uniprot_id, prop)].append(value)
    target_annotation = [
        {"target_key": accession, "source_db": SOURCE_DB, "property": prop,
         "ordinal": ordinal, "value": value}
        for (accession, prop), values in sorted(distinct.items())
        for ordinal, value in enumerate(values)]

    in_vivo = []
    for row in invivo[invivo.probe_ix.isin(keyed)].itertuples():
        doses = parse_dose(row.dose) or [(None, None, None)]
        for value, unit, route in doses:
            in_vivo.append({"inchikey": key_of[row.probe_ix], "organism": row.organism,
                            "dose_value": as_number(value), "dose_unit": unit or "",
                            "route": route or "", "dose_raw": row.dose,
                            "source_db": SOURCE_DB, "source": name_of[row.probe_ix]})

    references = []
    for row in reference[reference.probe_ix.isin(keyed)].itertuples():
        prefix, xref = reference_source(row.ref)
        references.append({"inchikey": key_of[row.probe_ix], "xref_id": prefix or "",
                           "source_xref": xref, "raw": row.ref})

    validations_per_probe = (validation.groupby("probe_ix").size()
                             if validation is not None else {})
    targets_per_probe = target.groupby("probe_ix").size()
    skipped = [{"name": row.name,
                "reason": "no InChIKey in the export and no structure resolved",
                "portal_path": path_of[row.probe_ix],
                "targets": targets_per_probe.get(row.probe_ix, 0),
                "validations": (validations_per_probe.get(row.probe_ix, 0)
                                if len(validations_per_probe) else "")}
               for row in probe[probe.inchikey == ""].itertuples()]

    return {"compound_annotation": annotation, "target_annotation": target_annotation,
            "in_vivo": in_vivo, "reference": references, "skipped_compound": skipped}


def write_table(path, rows, columns):
    """One list of dicts -> one TSV.

    csv.writer, never '\\t'.join: assayDesc carries tabs, newlines and quotes,
    which a naive join spreads across rows. Every value goes through as_text, so
    a NaN never lands as the string 'nan'.
    """
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    # a column with one missing value is float64 in pandas, which would turn
    # canSAR_ID 1354531 into '1354531.0'
    rows = rows.to_dict("records") if hasattr(rows, "to_dict") else list(rows)
    with open(path, "w", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t", quoting=csv.QUOTE_MINIMAL)
        writer.writerow(columns)
        for row in rows:
            writer.writerow([cell(row.get(column)) for column in columns])
    return len(rows)


def cell(value):
    """One value as it is written: a whole float loses its decimal point."""
    return as_number(value) if isinstance(value, float) else as_text(value)


def write_staging(out, tables, leftovers=None, quarantined=None, probes=None):
    out = Path(out)
    everything = dict(tables)
    everything.update(leftovers or {})
    if quarantined is not None:
        everything["quarantine"] = quarantined
    written = {}
    for name in STAGING_FILES + EXTRA_FILES:
        written[name] = write_table(out / f"{name}.tsv", everything.get(name, []),
                                    COLUMNS[name])
    summary = report(tables, leftovers, quarantined, probes=probes)
    summary["written"] = written
    (out / "report.json").write_text(json.dumps(summary, indent=2))
    return summary


def report(tables, leftovers=None, quarantined=None, probes=None):
    """Row counts per file, and the reconciliation that every input record
    reached exactly one output.

    Every probe is either a compound row or a skipped one. Every number found
    under a writable compound is either written to bioactivity.tsv or held back
    in quarantine.tsv, never both (D15).
    """
    leftovers = leftovers or {}
    quarantined = [] if quarantined is None else quarantined
    counted = dict(tables)
    counted.update(leftovers)
    counted["quarantine"] = quarantined
    written = {name: len(counted.get(name, [])) for name in STAGING_FILES + EXTRA_FILES}

    compounds, skipped = written["compound"], written["skipped_compound"]
    counts = {
        "probes in": len(probes) if probes is not None else compounds + skipped,
        "compounds written": compounds,
        "compounds skipped": skipped,
        "targets written": written["target"],
        "bioactivity rows written": written["bioactivity"],
        "rows held for curation": written["quarantine"],
        "rows with no unit": sum(1 for r in tables.get("bioactivity", [])
                                 if not r.get("unit")),
        "in vivo rows written": written["in_vivo"],
        "reference rows written": written["reference"],
        "unsuitable rows written": written["unsuitable"],
    }
    problems = []
    if probes is not None and compounds + skipped != len(probes):
        problems.append(f"{len(probes)} probes in, {compounds} written and "
                        f"{skipped} skipped")
    # the two files share these columns and no others, so this is the only
    # tuple that can be compared between them
    def measurement(row):
        return tuple(row.get(column) for column in
                     ("inchikey", "target_key", "value", "unit", "bioactivity_type"))

    held = {measurement(row) for row in quarantined}
    for row in tables.get("bioactivity", []):
        if measurement(row) in held:
            problems.append(f"held back and written at once: {measurement(row)}")
    return {"counts": counts, "written": written, "problems": problems}


# ------------------------------------------------------------------- entry

def preprocess(json_path, out, resolve=False):
    """The whole pipeline: export -> staging directory. Returns the report."""
    probes = load_export(json_path)
    if resolve:
        refresh_resolved_cache(probes)
    frames = flatten(probes)
    frames["probe"] = apply_resolved(frames["probe"],
                                    resolve_missing_structures(frames["probe"]))
    probe, target = frames["probe"], frames["target"]
    keyed = set(probe.loc[probe.inchikey != "", "probe_ix"])
    mine = target[target.probe_ix.isin(keyed)]

    rows, quarantined = split_bioactivity(probe, target, frames["validation"],
                                          frames["reference"])
    tables = {
        "compound": build_compound(probe, frames["chembl"]),
        "target": build_target(mine),
        "uniprot": build_uniprot(mine),
        "bioactivity": rows,
        "unsuitable": build_unsuitable(probe, frames["chembl"], frames["reference"]),
    }
    leftovers = build_leftovers(probe, target, frames["invivo"], frames["reference"],
                               frames["control"], frames["validation"])
    return write_staging(out, tables, leftovers, quarantined, probes=probes)


def refuse_existing(db_path, force):
    """create=True needs empty tables, so an existing file has to go first."""
    path = Path(db_path)
    if str(path) == ":memory:" or not path.exists():
        return
    if not force:
        raise FileExistsError(
            f"{path} already exists. create=True needs empty tables, so pass "
            f"--force to replace it or delete it yourself"
        )
    path.unlink()


def populate(db_path, out, force=False):
    """Build a fresh database from the staging directory, through loader/."""
    from loader import load
    from probedb import ProbeDB

    refuse_existing(db_path, force)
    db = ProbeDB(str(db_path), create=True)
    result = load(db, out, source=SOURCE_DB)
    return db, result


def main(argv=None):
    parser = argparse.ArgumentParser(
        description=__doc__.split("Decided:")[0].rstrip(),
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )
    parser.add_argument("--json", type=Path, required=True,
                        help="the portal export to read")
    parser.add_argument("--out", type=Path, required=True,
                        help="the staging directory to write, created if missing")
    parser.add_argument("--db", type=Path, help="also build a database from the result")
    parser.add_argument("--force", action="store_true",
                        help="replace the database file if it already exists")
    parser.add_argument("--resolve", action="store_true",
                        help="re-query PubChem for the probes with no structure "
                             "and rewrite resolved_structures.tsv")
    parser.add_argument("--validate", action="store_true",
                        help="run loader.validate on the result and print the problems. "
                             "Expect one 'no relation, kept anyway' per row without an "
                             "operator, about 3500, and none of them is a hard error")
    args = parser.parse_args(argv)
    if args.db:
        refuse_existing(args.db, args.force)      # before doing all the work
    check_vocabularies()

    summary = preprocess(args.json, args.out, resolve=args.resolve)
    logger.info("\n" + "\n".join(f"{v:>7}  {k}" for k, v in summary["counts"].items()))
    logger.info("\n" + "\n".join(f"{n:>7}  {name}.tsv"
                                 for name, n in summary["written"].items()))
    for problem in summary["problems"]:
        logger.error(problem)

    if args.validate:
        from loader import validate
        problems = validate(args.out)
        hard = [p for p in problems if not p.endswith(("ignored", "kept anyway"))]
        logger.info(f"validate: {len(problems)} problems, "
                    f"{len(hard)} of them hard errors")
        for problem in hard[:20]:
            logger.error(problem)

    if args.db:
        db, result = populate(args.db, args.out, force=args.force)
        logger.info(f"loaded into {args.db}\n"
                    + db.counts().to_string(index=False)
                    + f"\nduplicates skipped by the loader: "
                      f"{result['duplicates_skipped']}")
        db.close()
    return 0 if not summary["problems"] else 1


if __name__ == "__main__":
    sys.exit(main())
