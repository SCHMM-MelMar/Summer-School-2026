#!/usr/bin/env python3
"""
build_tables.py
================

Turns a SPARK-format supplementary table (e.g. pnas.1911532117.sd06.xlsx,
columns: IndexSpark1.01 | ChEMBLID | GeneSymbol | pActMean | ChEMBLPrefName)
into the four reference tables Evgenia is using:

    1. compounds.tsv          -> inchikey, smiles, chembl_id, name
    2. uniprot.tsv            -> uniprot_id, target_key, hgnc, species
    3. target.tsv  -> target_key, type, name
    4. bioactivity.tsv        -> inchikey, target_key, moa, bioactivity_type,
                                  relation, value, unit, assay_type,
                                  assay_description, cell_line, concentration,
                                  concentration_unit, source_db, source,
                                  source_xref, xref_id

Concept / logic (mirrors the example screenshots):
  * "uniprot" is a long table: ONE ROW PER (uniprot_id, target_key) PAIR.
    A single-protein target (e.g. PLK1) has exactly one row.
    A multi-subunit target (a "family" or "complex", e.g. "PARP 1, 2 and 3"
    or "CDK1/cyclin B1") has one row PER MEMBER PROTEIN, all sharing the
    same target_key.
  * "target" is the short lookup table: ONE ROW PER target_key,
    with its type (protein / family / complex) and display name.
  * "compounds" is one row per unique ChEMBL compound, enriched with its
    canonical SMILES and InChIKey.
  * "bioactivity" is one row per (compound, target) measurement, linking
    compounds.inchikey to targets.target_key. The SPARK sheet only gives
    us ChEMBLID + GeneSymbol + pActMean per row, so several bioactivity
    columns simply don't exist in the source data (moa, assay_type,
    assay_description, cell_line, concentration*) and are left blank /
    constant on purpose -- see "Manual overrides" below to fill them in
    if you have that information elsewhere.

Usage:
    python build_tables.py your_spark_file.xlsx --out ./out
    (add --sheet "SheetName" if your data isn't on the first sheet)

Outputs (in --out, default "./out"):
    compounds.csv
    uniprot.csv
    target.csv
    bioactivity.csv
    reference_tables.xlsx   (all four as separate sheets)
    api_cache.json          (cache of every API lookup -> safe to re-run)
"""

import argparse
import json
import sys
import time
from pathlib import Path

import pandas as pd
import requests

CHEMBL_MOLECULE_URL = "https://www.ebi.ac.uk/chembl/api/data/molecule.json"
CHEMBL_MOLECULE_SINGLE_URL = "https://www.ebi.ac.uk/chembl/api/data/molecule/{}.json"
CHEMBL_ONLY_FIELDS = "molecule_chembl_id,pref_name,molecule_structures"
UNIPROT_SEARCH_URL = "https://rest.uniprot.org/uniprotkb/search"

CHUNK_SIZE_CHEMBL = 20      # molecule_chembl_id__in chunk size (smaller
                            # chunks + the per-ID fallback below make the
                            # flaky bulk endpoint reliable)
CHUNK_SIZE_UNIPROT = 40     # gene symbols per UniProt "OR" query
REQUEST_TIMEOUT = 30
RETRIES = 3
SLEEP_BETWEEN_RETRIES = 2
SLEEP_BETWEEN_REQUESTS = 0.34  # ~3 requests/sec, polite to the public APIs

# Preferred organism order when several species match the same gene symbol.
ORGANISM_PRIORITY = {
    "Homo sapiens": 0,
    "Mus musculus": 1,
    "Rattus norvegicus": 2,
}

# ---------------------------------------------------------------------------
# Manual overrides
# ---------------------------------------------------------------------------
# Deciding "family" vs "complex" for multi-gene targets is a biological call,
# not something derivable from the spreadsheet alone. By default every
# multi-gene target_key is classified as "complex". List target_keys here
# that should instead be a "family" (e.g. isoenzyme groups like PARP1/2/3,
# adrenergic receptor subtypes, carbonic anhydrases, etc.).
FAMILY_OVERRIDES = {
    # "ADRA1A/ADRA1B/ADRA1D": "family",
}

# Optional friendly display names / target_key overrides, keyed by the
# semicolon-joined GeneSymbol string exactly as it appears in the source
# spreadsheet. Anything not listed here falls back to an automatic
# slash-joined target_key (e.g. "FKBP1A;MTOR" -> "FKBP1A/MTOR").
TARGET_KEY_OVERRIDES = {
    # "CCNB1;CDK1": "CDK1/cyclin B1",
}

# ---------------------------------------------------------------------------
# Bioactivity table constants
# ---------------------------------------------------------------------------
# The SPARK sheet reports one pre-aggregated "pActMean" per (compound,
# target) pair -- the mean of -log10(molar activity) across whatever mix of
# IC50/Ki/EC50-style assays ChEMBL had for that pair. Because it's a mean
# across potentially different assay types, we deliberately do NOT
# back-convert it into a single "IC50 in nM" value (that would overstate
# precision). It's reported as-is, in its native pAct unit.
BIOACTIVITY_TYPE = "pActMean"
BIOACTIVITY_UNIT = "-log10(M)"
BIOACTIVITY_RELATION = None  # it's a mean, not a single censored measurement
SOURCE_DB = "SPARK-PNAS"
# Fill in the exact citation/DOI you want recorded against every row, e.g.
# "Cortes-Ciriano et al., PNAS 2020, 10.1073/pnas.1911532117".
SOURCE_CITATION = None


def load_cache(cache_path: Path) -> dict:
    if cache_path.exists():
        with open(cache_path, "r") as f:
            return json.load(f)
    return {"chembl": {}, "uniprot": {}}


def save_cache(cache_path: Path, cache: dict) -> None:
    with open(cache_path, "w") as f:
        json.dump(cache, f, indent=2, sort_keys=True)


def chunked(seq, size):
    for i in range(0, len(seq), size):
        yield seq[i : i + size]


def fetch_chembl_batch(chembl_ids, cache):
    """Populate cache['chembl'][chembl_id] = {smiles, inchikey, pref_name}.

    Tries the bulk molecule.json?molecule_chembl_id__in=... lookup first
    (fast). ChEMBL's public API intermittently 500s on that endpoint under
    load, so anything not resolved by the bulk call is retried one ID at a
    time against /molecule/{id}.json, which is more reliable. Both calls
    request only the 3 fields we need (?only=...) -- some individual
    ChEMBL records otherwise 500 permanently because an unrelated field on
    that record crashes the default full-record serializer.
    """
    todo = [c for c in chembl_ids if c not in cache["chembl"]]
    for chunk in chunked(todo, CHUNK_SIZE_CHEMBL):
        params = {
            "molecule_chembl_id__in": ",".join(chunk),
            "only": CHEMBL_ONLY_FIELDS,
            "format": "json",
            "limit": len(chunk),
        }
        data = _get_with_retry(CHEMBL_MOLECULE_URL, params, quiet_errors=True)
        found = set()
        if data:
            for mol in data.get("molecules", []):
                cid = mol.get("molecule_chembl_id")
                struct = mol.get("molecule_structures") or {}
                cache["chembl"][cid] = {
                    "smiles": struct.get("canonical_smiles"),
                    "inchikey": struct.get("standard_inchi_key"),
                    "pref_name": mol.get("pref_name"),
                }
                found.add(cid)
        time.sleep(SLEEP_BETWEEN_REQUESTS)

        missing = [c for c in chunk if c not in found]
        if missing:
            if not data:
                print(f"  [info] bulk lookup failed for this chunk of {len(chunk)}; "
                      f"falling back to one-by-one lookups (slower, more reliable)")
            for cid in missing:
                mol = _get_with_retry(
                    CHEMBL_MOLECULE_SINGLE_URL.format(cid),
                    {"format": "json", "only": CHEMBL_ONLY_FIELDS},
                )
                struct = (mol or {}).get("molecule_structures") or {}
                cache["chembl"][cid] = {
                    "smiles": struct.get("canonical_smiles"),
                    "inchikey": struct.get("standard_inchi_key"),
                    "pref_name": (mol or {}).get("pref_name"),
                }
                time.sleep(SLEEP_BETWEEN_REQUESTS)

        resolved_now = sum(1 for c in chunk if cache["chembl"].get(c, {}).get("smiles"))
        total_resolved = sum(1 for v in cache["chembl"].values() if v.get("smiles"))
        print(f"  ChEMBL: resolved {resolved_now}/{len(chunk)} in this chunk "
              f"({total_resolved} total so far)")


def fetch_uniprot_batch(gene_symbols, cache):
    """Populate cache['uniprot'][gene_symbol] = list of candidate hits."""
    todo = [g for g in gene_symbols if g not in cache["uniprot"]]
    for chunk in chunked(todo, CHUNK_SIZE_UNIPROT):
        gene_query = " OR ".join(f'gene:{g}' for g in chunk)
        query = f"({gene_query}) AND (organism_id:9606 OR organism_id:10090 OR organism_id:10116)"
        params = {
            "query": query,
            "fields": "accession,gene_names,organism_name,protein_name,reviewed",
            "format": "json",
            "size": 500,
        }
        data = _get_with_retry(UNIPROT_SEARCH_URL, params)
        by_gene = {g: [] for g in chunk}
        if data:
            for entry in data.get("results", []):
                genes = entry.get("genes") or []
                gene_names = set()
                for g in genes:
                    if g.get("geneName", {}).get("value"):
                        gene_names.add(g["geneName"]["value"])
                    for syn in g.get("synonyms", []) or []:
                        if syn.get("value"):
                            gene_names.add(syn["value"])
                organism = (entry.get("organism") or {}).get("scientificName")
                accession = entry.get("primaryAccession")
                reviewed = entry.get("entryType", "").startswith("UniProtKB reviewed")
                full_name = (
                    (entry.get("proteinDescription") or {})
                    .get("recommendedName", {})
                    .get("fullName", {})
                    .get("value")
                )
                for wanted in chunk:
                    if wanted in gene_names or wanted.upper() in {n.upper() for n in gene_names}:
                        by_gene[wanted].append(
                            {
                                "uniprot_id": accession,
                                "organism": organism,
                                "reviewed": reviewed,
                                "full_name": full_name,
                            }
                        )
        for g in chunk:
            cache["uniprot"][g] = by_gene[g]
        time.sleep(SLEEP_BETWEEN_REQUESTS)
        print(f"  UniProt: queried {len(chunk)} gene symbols "
              f"({sum(1 for v in cache['uniprot'].values() if v)} resolved so far)")


def _get_with_retry(url, params, quiet_errors=False):
    """quiet_errors=True skips the per-attempt/give-up prints -- used for the
    bulk ChEMBL lookup, where a failure is expected to sometimes happen and
    is silently handled by the per-ID fallback in fetch_chembl_batch, so
    logging it would just be noise."""
    for attempt in range(1, RETRIES + 1):
        try:
            resp = requests.get(url, params=params, timeout=REQUEST_TIMEOUT)
            if resp.status_code == 200:
                return resp.json()
            if not quiet_errors:
                print(f"  [warn] HTTP {resp.status_code} for {url} (attempt {attempt}/{RETRIES})")
        except requests.RequestException as e:
            if not quiet_errors:
                print(f"  [warn] {e} (attempt {attempt}/{RETRIES})")
        time.sleep(SLEEP_BETWEEN_RETRIES * attempt)
    if not quiet_errors:
        print(f"  [error] giving up on {url} with params={params}")
    return None


def pick_best_uniprot_hit(hits):
    """Prefer reviewed (Swiss-Prot) entries, then Homo sapiens > Mus > Rattus."""
    if not hits:
        return None

    def sort_key(h):
        organism_rank = ORGANISM_PRIORITY.get(h["organism"], 99)
        reviewed_rank = 0 if h["reviewed"] else 1
        return (reviewed_rank, organism_rank)

    return sorted(hits, key=sort_key)[0]


def make_target_key(gene_symbol_field: str) -> str:
    if gene_symbol_field in TARGET_KEY_OVERRIDES:
        return TARGET_KEY_OVERRIDES[gene_symbol_field]
    genes = [g.strip() for g in gene_symbol_field.split(";") if g.strip()]
    # de-duplicate while preserving order (source data has literal dupes,
    # e.g. "APOBEC3A;APOBEC3A")
    seen = []
    for g in genes:
        if g not in seen:
            seen.append(g)
    return "/".join(seen)


def build_tables(df: pd.DataFrame, cache: dict):
    df = df.dropna(subset=["ChEMBLID", "GeneSymbol"]).copy()
    df["GeneSymbol"] = df["GeneSymbol"].astype(str).str.strip()
    df["ChEMBLID"] = df["ChEMBLID"].astype(str).str.strip()

    # ---- compounds -------------------------------------------------------
    compound_rows = (
        df[["ChEMBLID", "ChEMBLPrefName"]]
        .drop_duplicates(subset=["ChEMBLID"])
        .rename(columns={"ChEMBLID": "chembl_id", "ChEMBLPrefName": "sheet_name"})
    )
    chembl_ids = compound_rows["chembl_id"].tolist()
    print(f"Fetching ChEMBL data for {len(chembl_ids)} unique compounds...")
    fetch_chembl_batch(chembl_ids, cache)

    def compound_name(row):
        info = cache["chembl"].get(row["chembl_id"], {})
        return info.get("pref_name") or row["sheet_name"]

    compound_rows["name"] = compound_rows.apply(compound_name, axis=1)
    compound_rows["smiles"] = compound_rows["chembl_id"].map(
        lambda c: cache["chembl"].get(c, {}).get("smiles")
    )
    compound_rows["inchikey"] = compound_rows["chembl_id"].map(
        lambda c: cache["chembl"].get(c, {}).get("inchikey")
    )
    compounds = compound_rows[["inchikey", "smiles", "chembl_id", "name"]].sort_values("chembl_id")

    # ---- targets / target_dictionary -------------------------------------
    unique_gene_fields = sorted(df["GeneSymbol"].unique())
    target_key_map = {gf: make_target_key(gf) for gf in unique_gene_fields}

    all_genes = sorted({g for tk in target_key_map.values() for g in tk.split("/")})
    print(f"Fetching UniProt data for {len(all_genes)} unique gene symbols...")
    fetch_uniprot_batch(all_genes, cache)

    target_rows = []
    dict_rows = []
    seen_target_keys = set()

    for gene_field, target_key in target_key_map.items():
        member_genes = target_key.split("/")
        is_multi = len(member_genes) > 1

        member_infos = []
        for gene in member_genes:
            best = pick_best_uniprot_hit(cache["uniprot"].get(gene, []))
            member_infos.append((gene, best))
            target_rows.append(
                {
                    "uniprot_id": best["uniprot_id"] if best else None,
                    "target_key": target_key,
                    "hgnc": gene,
                    "species": best["organism"] if best else None,
                }
            )

        if target_key not in seen_target_keys:
            seen_target_keys.add(target_key)
            if is_multi:
                ttype = FAMILY_OVERRIDES.get(target_key, "complex")
                names = [info["full_name"] for _, info in member_infos if info and info["full_name"]]
                name = " / ".join(names) if names else target_key
            else:
                ttype = "protein"
                gene, info = member_infos[0]
                name = info["full_name"] if info and info["full_name"] else gene
            dict_rows.append({"target_key": target_key, "type": ttype, "name": name})

    targets = pd.DataFrame(target_rows).drop_duplicates().sort_values(["target_key", "hgnc"])
    target_dictionary = pd.DataFrame(dict_rows).sort_values("target_key")

    # ---- bioactivity -------------------------------------------------------
    bioactivity = build_bioactivity(df, compounds, target_key_map)

    return compounds, targets, target_dictionary, bioactivity


def build_bioactivity(df: pd.DataFrame, compounds: pd.DataFrame, target_key_map: dict) -> pd.DataFrame:
    """One row per (compound, target) measurement from the SPARK sheet.

    inchikey/target_key are the real join keys back to the compounds and
    targets tables; everything the source sheet doesn't tell us
    (moa, assay_type, assay_description, cell_line, concentration*) is left
    as None on purpose rather than guessed.
    """
    inchikey_by_chembl = dict(zip(compounds["chembl_id"], compounds["inchikey"]))

    rows = []
    for _, row in df.iterrows():
        chembl_id = row["ChEMBLID"]
        gene_field = row["GeneSymbol"]
        rows.append(
            {
                "inchikey": inchikey_by_chembl.get(chembl_id),
                "target_key": target_key_map.get(gene_field),
                "moa": None,  # not present in the SPARK sheet
                "bioactivity_type": BIOACTIVITY_TYPE,
                "relation": BIOACTIVITY_RELATION,
                "value": row.get("pActMean"),
                "unit": BIOACTIVITY_UNIT,
                "assay_type": None,  # not present in the SPARK sheet
                "assay_description": None,  # not present in the SPARK sheet
                "cell_line": None,  # not present in the SPARK sheet
                "concentration": None,  # not present in the SPARK sheet
                "concentration_unit": None,  # not present in the SPARK sheet
                "source_db": SOURCE_DB,
                "source": SOURCE_CITATION,
                "source_xref": f"IndexSpark1.01:{row['IndexSpark1.01']}",
                "xref_id": None,
            }
        )

    columns = [
        "inchikey", "target_key", "moa", "bioactivity_type", "relation", "value",
        "unit", "assay_type", "assay_description", "cell_line", "concentration",
        "concentration_unit", "source_db", "source", "source_xref", "xref_id",
    ]
    return pd.DataFrame(rows, columns=columns)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("input_xlsx", help="Path to your SPARK-format .xlsx")
    parser.add_argument(
        "--sheet", default=None,
        help="Sheet name to read. Default: the first sheet in the file.",
    )
    parser.add_argument("--out", default="./out", help="Output directory")
    args = parser.parse_args()

    out_dir = Path(args.out)
    out_dir.mkdir(parents=True, exist_ok=True)
    cache_path = out_dir / "api_cache.json"
    cache = load_cache(cache_path)

    sheet = args.sheet if args.sheet is not None else 0  # 0 = first sheet
    print(f"Reading {args.input_xlsx} [sheet={sheet!r}] ...")
    df = pd.read_excel(args.input_xlsx, sheet_name=sheet)

    required = {"ChEMBLID", "GeneSymbol", "pActMean", "ChEMBLPrefName", "IndexSpark1.01"}
    missing = required - set(df.columns)
    if missing:
        sys.exit(
            f"Missing expected column(s) {sorted(missing)} in sheet {sheet!r}. "
            f"Columns found: {list(df.columns)}. "
            f"Pass --sheet <name> if your data isn't on the first sheet."
        )

    try:
        compounds, targets, target_dictionary, bioactivity = build_tables(df, cache)
    finally:
        save_cache(cache_path, cache)  # keep whatever we fetched even on error/interrupt

    compounds.to_csv(out_dir / "compounds.tsv", index=False)
    targets.to_csv(out_dir / "uniprot.tsv", index=False)
    target_dictionary.to_csv(out_dir / "targets.tsv", index=False)
    bioactivity.to_csv(out_dir / "bioactivity.tsv", index=False)

    with pd.ExcelWriter(out_dir / "reference_tables.xlsx") as writer:
        compounds.to_excel(writer, sheet_name="compounds", index=False)
        targets.to_excel(writer, sheet_name="targets", index=False)
        target_dictionary.to_excel(writer, sheet_name="target_dictionary", index=False)
        bioactivity.to_excel(writer, sheet_name="bioactivity", index=False)

    print("\nDone.")
    print(f"  compounds:          {len(compounds)} rows")
    print(f"  targets:            {len(targets)} rows")
    print(f"  target_dictionary:  {len(target_dictionary)} rows")
    print(f"  bioactivity:        {len(bioactivity)} rows")
    missing_smiles = compounds["smiles"].isna().sum()
    missing_uniprot = targets["uniprot_id"].isna().sum()
    if missing_smiles == len(compounds) or missing_uniprot == len(targets):
        print(
            "\n  [!] ALL compound/target lookups came back empty -- this almost always means "
            "www.ebi.ac.uk and/or rest.uniprot.org were unreachable from this network "
            "(corporate/institutional proxy blocking them), not a bug. See the note at the "
            "top of this script for how to work around it."
        )
    elif missing_smiles or missing_uniprot:
        if missing_smiles:
            print(f"  [!] {missing_smiles} compounds have no SMILES (ChEMBL id not found/withdrawn "
                  f"even after the per-ID fallback -- see compounds.csv for which ones)")
        if missing_uniprot:
            print(f"  [!] {missing_uniprot} target rows have no uniprot_id (gene symbol not matched)")
    print(f"\nOutputs written to: {out_dir.resolve()}")


if __name__ == "__main__":
    main()