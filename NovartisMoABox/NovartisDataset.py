from __future__ import annotations

import argparse
import csv
import re
import time
from pathlib import Path
from typing import Iterable

import pandas as pd
import requests


DEFAULT_EXCEL = "March2020_NIBRmoabox-data_OAK.xlsx"
DEFAULT_SHEET = "Report"
DEFAULT_OUTPUT_DIR = Path(__file__).resolve().parent / "custom_data"
DEFAULT_CHEMBL_OUTPUT = DEFAULT_OUTPUT_DIR / "bioactivity.tsv"
SOURCE_DB = "NIBR MOA Box"
SPECIES = "Homo sapiens"
CHEMBL_BASE_URL = "https://www.ebi.ac.uk/chembl/api/data"


BIOACTIVITY_COLUMNS = [
    "inchikey",
    "target_key",
    "moa",
    "bioactivity_type",
    "relation",
    "value",
    "unit",
    "assay_type",
    "assay_description",
    "cell_line",
    "concentration",
    "concentration_unit",
    "source_db",
    "source",
    "source_xref",
    "xref_id",
]

COMPOUND_COLUMNS = ["inchikey", "smiles", "chembl_id", "name"]
TARGET_COLUMNS = ["target_key", "type", "name"]
UNIPROT_COLUMNS = ["uniprot_id", "target_key", "hgnc", "species"]


def clean_value(value: object) -> str:
    """Return a TSV-friendly empty string for missing values."""
    if pd.isna(value):
        return ""
    if isinstance(value, float) and value.is_integer():
        return str(int(value))
    return str(value).strip()


def first_present(row: pd.Series, columns: Iterable[str]) -> str:
    for column in columns:
        if column in row and pd.notna(row[column]):
            return clean_value(row[column])
    return ""


def first_present_dict(record: dict, keys: Iterable[str]) -> str:
    for key in keys:
        value = record.get(key)
        if value not in (None, ""):
            return clean_value(value)
    return ""


def indexed_columns(df: pd.DataFrame, prefix: str) -> dict[int, str]:
    pattern = re.compile(rf"^{re.escape(prefix)}\[(\d+)\]$")
    found: dict[int, str] = {}
    for column in df.columns:
        match = pattern.match(column)
        if match:
            found[int(match.group(1))] = column
    return found


def load_report(excel_path: Path, sheet_name: str) -> pd.DataFrame:
    return pd.read_excel(excel_path, sheet_name=sheet_name)


def make_compound(df: pd.DataFrame) -> pd.DataFrame:
    chembl_columns = [column for _, column in sorted(indexed_columns(df, "chembl_ids").items())]
    records = []

    for _, row in df.iterrows():
        records.append(
            {
                "inchikey": clean_value(row.get("inchi_key")),
                "smiles": clean_value(row.get("smiles")),
                "chembl_id": first_present(row, chembl_columns),
                "name": "",
            }
        )

    return (
        pd.DataFrame(records, columns=COMPOUND_COLUMNS)
        .query("inchikey != ''")
        .drop_duplicates(subset=["inchikey"], keep="first")
        .sort_values("inchikey")
    )


def make_target(df: pd.DataFrame) -> pd.DataFrame:
    gene_columns = [column for _, column in sorted(indexed_columns(df, "gene_symbols").items())]
    target_keys = sorted(
        {
            clean_value(value)
            for column in gene_columns
            for value in df[column].dropna().tolist()
            if clean_value(value)
        }
    )

    records = [
        {"target_key": target_key, "type": "protein", "name": target_key}
        for target_key in target_keys
    ]
    return pd.DataFrame(records, columns=TARGET_COLUMNS)


def make_uniprot(df: pd.DataFrame) -> pd.DataFrame:
    gene_columns = [column for _, column in sorted(indexed_columns(df, "gene_symbols").items())]
    uniprot_columns = [
        column
        for column in df.columns
        if "uniprot" in column.lower() or "swiss" in column.lower()
    ]
    records_by_target: dict[str, dict[str, str]] = {}

    for _, row in df.iterrows():
        for gene_column in gene_columns:
            target_key = clean_value(row.get(gene_column))
            if not target_key or target_key in records_by_target:
                continue

            records_by_target[target_key] = {
                "uniprot_id": first_present(row, uniprot_columns),
                "target_key": target_key,
                "hgnc": target_key,
                "species": SPECIES,
            }

    return pd.DataFrame(
        [records_by_target[key] for key in sorted(records_by_target)],
        columns=UNIPROT_COLUMNS,
    )


def make_bioactivity(df: pd.DataFrame, source_name: str) -> pd.DataFrame:
    gene_columns = indexed_columns(df, "gene_symbols")
    measure_columns = {
        "target_score": indexed_columns(df, "target_scores"),
        "selectivity": indexed_columns(df, "selectivities"),
        "active_strength": indexed_columns(df, "active_strengths"),
    }
    records = []

    for _, row in df.iterrows():
        inchikey = clean_value(row.get("inchi_key"))
        if not inchikey:
            continue

        for index, gene_column in sorted(gene_columns.items()):
            target_key = clean_value(row.get(gene_column))
            if not target_key:
                continue

            for measure_name, columns_by_index in measure_columns.items():
                measure_column = columns_by_index.get(index)
                if not measure_column or pd.isna(row.get(measure_column)):
                    continue

                records.append(
                    {
                        "inchikey": inchikey,
                        "target_key": target_key,
                        "moa": clean_value(row.get("moa")),
                        "bioactivity_type": measure_name,
                        "relation": "=",
                        "value": clean_value(row.get(measure_column)),
                        "unit": "score",
                        "assay_type": "",
                        "assay_description": "",
                        "cell_line": "",
                        "concentration": "",
                        "concentration_unit": "",
                        "source_db": SOURCE_DB,
                        "source": source_name,
                        "source_xref": clean_value(row.get("moabox_id")),
                        "xref_id": "",
                    }
                )

    return pd.DataFrame(records, columns=BIOACTIVITY_COLUMNS)


def all_chembl_ids(row: pd.Series, chembl_columns: Iterable[str]) -> list[str]:
    ids = []
    for column in chembl_columns:
        value = clean_value(row.get(column))
        if value and value not in ids:
            ids.append(value)
    return ids


def make_workbook_compound_maps(df: pd.DataFrame) -> tuple[dict[str, str], dict[str, str]]:
    chembl_columns = [column for _, column in sorted(indexed_columns(df, "chembl_ids").items())]
    chembl_to_inchikey: dict[str, str] = {}
    inchikey_to_moa: dict[str, str] = {}

    for _, row in df.iterrows():
        inchikey = clean_value(row.get("inchi_key"))
        if not inchikey:
            continue

        inchikey_to_moa.setdefault(inchikey, clean_value(row.get("moa")))
        for chembl_id in all_chembl_ids(row, chembl_columns):
            chembl_to_inchikey.setdefault(chembl_id, inchikey)

    return chembl_to_inchikey, inchikey_to_moa


def chembl_get(endpoint: str, params: dict | None = None, retries: int = 3) -> dict:
    url = f"{CHEMBL_BASE_URL}/{endpoint.lstrip('/')}.json"
    params = params or {}
    last_error: Exception | None = None

    for attempt in range(1, retries + 1):
        try:
            response = requests.get(url, params=params, timeout=60)
            response.raise_for_status()
            return response.json()
        except (requests.RequestException, ValueError) as error:
            last_error = error
            if attempt == retries:
                break
            time.sleep(2 * attempt)

    raise RuntimeError(f"ChEMBL request failed for {url}: {last_error}") from last_error


def iter_chembl_records(endpoint: str, collection_key: str, params: dict) -> Iterable[dict]:
    offset = 0
    limit = int(params.get("limit", 1000))

    while True:
        page_params = {**params, "limit": limit, "offset": offset}
        payload = chembl_get(endpoint, page_params)
        records = payload.get(collection_key, [])
        if not records:
            break

        yield from records

        page_meta = payload.get("page_meta", {})
        if not page_meta.get("next") or len(records) < limit:
            break
        offset += limit


def lookup_chembl_ids_by_inchikey(inchikeys: Iterable[str]) -> dict[str, str]:
    inchikey_to_chembl: dict[str, str] = {}

    for inchikey in inchikeys:
        payload = chembl_get(
            "molecule",
            {"molecule_structures__standard_inchi_key": inchikey, "limit": 1},
        )
        molecules = payload.get("molecules", [])
        if molecules:
            chembl_id = clean_value(molecules[0].get("molecule_chembl_id"))
            if chembl_id:
                inchikey_to_chembl[inchikey] = chembl_id

    return inchikey_to_chembl


def extract_target_key(target_record: dict, fallback: str) -> str:
    for component in target_record.get("target_components", []) or []:
        for synonym in component.get("target_component_synonyms", []) or []:
            if synonym.get("syn_type") == "GENE_SYMBOL":
                symbol = clean_value(synonym.get("component_synonym"))
                if symbol:
                    return symbol

        description = clean_value(component.get("component_description"))
        if description:
            return description

    return first_present_dict(target_record, ["pref_name", "target_pref_name", "target_chembl_id"]) or fallback


def make_target_key_resolver(enable_lookup: bool = True) -> callable:
    cache: dict[str, str] = {}

    def resolve(target_chembl_id: str, fallback: str = "") -> str:
        if not enable_lookup:
            return fallback or target_chembl_id
        if not target_chembl_id:
            return fallback
        if target_chembl_id not in cache:
            try:
                payload = chembl_get("target", {"target_chembl_id": target_chembl_id, "limit": 1})
                targets = payload.get("targets", [])
                cache[target_chembl_id] = (
                    extract_target_key(targets[0], target_chembl_id) if targets else target_chembl_id
                )
            except RuntimeError:
                cache[target_chembl_id] = fallback or target_chembl_id
        return cache[target_chembl_id]

    return resolve


def make_chembl_bioactivity(
    df: pd.DataFrame,
    max_compounds: int | None = None,
    max_activities: int | None = None,
    fallback_to_inchikey_lookup: bool = True,
    resolve_targets: bool = True,
) -> pd.DataFrame:
    chembl_to_inchikey, inchikey_to_moa = make_workbook_compound_maps(df)

    if fallback_to_inchikey_lookup:
        missing_inchikeys = [
            clean_value(value)
            for value in df["inchi_key"].dropna().tolist()
            if clean_value(value) and clean_value(value) not in set(chembl_to_inchikey.values())
        ]
        for inchikey, chembl_id in lookup_chembl_ids_by_inchikey(missing_inchikeys).items():
            chembl_to_inchikey.setdefault(chembl_id, inchikey)

    compound_ids = sorted(chembl_to_inchikey)
    if max_compounds is not None:
        compound_ids = compound_ids[:max_compounds]

    resolve_target_key = make_target_key_resolver(enable_lookup=resolve_targets)
    records = []

    for compound_number, compound_id in enumerate(compound_ids, start=1):
        if compound_number % 100 == 0:
            print(f"Queried {compound_number}/{len(compound_ids)} compounds...", flush=True)

        try:
            activities = iter_chembl_records(
                "activity",
                "activities",
                {"molecule_chembl_id": compound_id, "limit": 1000},
            )
            for activity in activities:
                inchikey = chembl_to_inchikey.get(compound_id, "")
                target_chembl_id = clean_value(activity.get("target_chembl_id"))
                records.append(
                    {
                        "inchikey": inchikey,
                        "target_key": resolve_target_key(
                            target_chembl_id,
                            first_present_dict(activity, ["target_pref_name", "target_chembl_id"]),
                        ),
                        "moa": inchikey_to_moa.get(inchikey, ""),
                        "bioactivity_type": clean_value(activity.get("standard_type")),
                        "relation": clean_value(activity.get("standard_relation")) or "=",
                        "value": clean_value(activity.get("standard_value")),
                        "unit": clean_value(activity.get("standard_units")) or "unspecified",
                        "assay_type": clean_value(activity.get("assay_type")),
                        "assay_description": clean_value(activity.get("assay_description")),
                        "cell_line": first_present_dict(
                            activity,
                            ["cell_chembl_id", "cell_id", "cell_line", "cell_type"],
                        ),
                        "concentration": "",
                        "concentration_unit": "",
                        "source_db": "ChEMBL",
                        "source": first_present_dict(activity, ["document_chembl_id", "src_id"]),
                        "source_xref": clean_value(activity.get("assay_chembl_id")),
                        "xref_id": clean_value(activity.get("activity_id")),
                    }
                )

                if max_activities is not None and len(records) >= max_activities:
                    return pd.DataFrame(records, columns=BIOACTIVITY_COLUMNS)
        except RuntimeError as error:
            print(f"Skipping {compound_id}: {error}", flush=True)
            continue

    return pd.DataFrame(records, columns=BIOACTIVITY_COLUMNS)


def iter_chembl_bioactivity_records(
    df: pd.DataFrame,
    max_compounds: int | None = None,
    max_activities: int | None = None,
    fallback_to_inchikey_lookup: bool = True,
    resolve_targets: bool = True,
) -> Iterable[dict[str, str]]:
    chembl_to_inchikey, inchikey_to_moa = make_workbook_compound_maps(df)

    if fallback_to_inchikey_lookup:
        missing_inchikeys = [
            clean_value(value)
            for value in df["inchi_key"].dropna().tolist()
            if clean_value(value) and clean_value(value) not in set(chembl_to_inchikey.values())
        ]
        for inchikey, chembl_id in lookup_chembl_ids_by_inchikey(missing_inchikeys).items():
            chembl_to_inchikey.setdefault(chembl_id, inchikey)

    compound_ids = sorted(chembl_to_inchikey)
    if max_compounds is not None:
        compound_ids = compound_ids[:max_compounds]

    emitted = 0
    resolve_target_key = make_target_key_resolver(enable_lookup=resolve_targets)

    for compound_number, compound_id in enumerate(compound_ids, start=1):
        if compound_number % 100 == 0:
            print(f"Queried {compound_number}/{len(compound_ids)} compounds...", flush=True)

        try:
            activities = iter_chembl_records(
                "activity",
                "activities",
                {"molecule_chembl_id": compound_id, "limit": 1000},
            )
            for activity in activities:
                inchikey = chembl_to_inchikey.get(compound_id, "")
                target_chembl_id = clean_value(activity.get("target_chembl_id"))
                yield {
                    "inchikey": inchikey,
                    "target_key": resolve_target_key(
                        target_chembl_id,
                        first_present_dict(activity, ["target_pref_name", "target_chembl_id"]),
                    ),
                    "moa": inchikey_to_moa.get(inchikey, ""),
                    "bioactivity_type": clean_value(activity.get("standard_type")),
                    "relation": clean_value(activity.get("standard_relation")) or "=",
                    "value": clean_value(activity.get("standard_value")),
                    "unit": clean_value(activity.get("standard_units")) or "unspecified",
                    "assay_type": clean_value(activity.get("assay_type")),
                    "assay_description": clean_value(activity.get("assay_description")),
                    "cell_line": first_present_dict(
                        activity,
                        ["cell_chembl_id", "cell_id", "cell_line", "cell_type"],
                    ),
                    "concentration": "",
                    "concentration_unit": "",
                    "source_db": "ChEMBL",
                    "source": first_present_dict(activity, ["document_chembl_id", "src_id"]),
                    "source_xref": clean_value(activity.get("assay_chembl_id")),
                    "xref_id": clean_value(activity.get("activity_id")),
                }
                emitted += 1

                if max_activities is not None and emitted >= max_activities:
                    return
        except RuntimeError as error:
            print(f"Skipping {compound_id}: {error}", flush=True)
            continue


def write_tsv(df: pd.DataFrame, output_path: Path) -> None:
    df.to_csv(output_path, sep="\t", index=False, lineterminator="\n")


def convert_excel_to_tsv(excel_path: Path, output_dir: Path, sheet_name: str) -> dict[str, int]:
    df = load_report(excel_path, sheet_name)
    output_dir.mkdir(parents=True, exist_ok=True)

    tables = {
        "bioactivity.tsv": make_bioactivity(df, excel_path.name),
        "compound.tsv": make_compound(df),
        "target.tsv": make_target(df),
        "uniprot.tsv": make_uniprot(df),
    }

    for filename, table in tables.items():
        write_tsv(table, output_dir / filename)

    return {filename: len(table) for filename, table in tables.items()}


def create_chembl_bioactivity_tsv(
    excel_path: Path,
    output_path: Path,
    sheet_name: str,
    max_compounds: int | None = None,
    max_activities: int | None = None,
    fallback_to_inchikey_lookup: bool = True,
    resolve_targets: bool = True,
) -> int:
    df = load_report(excel_path, sheet_name)
    output_path.parent.mkdir(parents=True, exist_ok=True)
    row_count = 0

    with output_path.open("w", newline="", encoding="utf-8") as output_file:
        writer = csv.DictWriter(
            output_file,
            fieldnames=BIOACTIVITY_COLUMNS,
            delimiter="\t",
            lineterminator="\n",
        )
        writer.writeheader()
        for record in iter_chembl_bioactivity_records(
            df,
            max_compounds=max_compounds,
            max_activities=max_activities,
            fallback_to_inchikey_lookup=fallback_to_inchikey_lookup,
            resolve_targets=resolve_targets,
        ):
            writer.writerow(record)
            row_count += 1
            if row_count % 1000 == 0:
                output_file.flush()
                print(f"Wrote {row_count} ChEMBL activity rows...", flush=True)

    return row_count


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Convert the Novartis MOA Box Excel workbook into normalized TSV files."
    )
    parser.add_argument(
        "excel",
        nargs="?",
        default=DEFAULT_EXCEL,
        help=f"Input Excel workbook. Defaults to {DEFAULT_EXCEL}.",
    )
    parser.add_argument(
        "-o",
        "--output-dir",
        default=str(DEFAULT_OUTPUT_DIR),
        help=f"Folder where the TSV files should be written. Defaults to {DEFAULT_OUTPUT_DIR}.",
    )
    parser.add_argument(
        "--sheet",
        default=DEFAULT_SHEET,
        help=f"Workbook sheet to read. Defaults to {DEFAULT_SHEET}.",
    )
    parser.add_argument(
        "--chembl-bioactivity",
        action="store_true",
        help="Write ChEMBL-derived bioactivity_chembl.tsv instead of the Excel-derived TSV set.",
    )
    parser.add_argument(
        "--chembl-output",
        default=str(DEFAULT_CHEMBL_OUTPUT),
        help=f"Output path for --chembl-bioactivity. Defaults to {DEFAULT_CHEMBL_OUTPUT}.",
    )
    parser.add_argument(
        "--max-compounds",
        type=int,
        default=None,
        help="Optional cap for ChEMBL compound queries, useful for testing.",
    )
    parser.add_argument(
        "--max-activities",
        type=int,
        default=None,
        help="Optional cap for ChEMBL activity rows, useful for testing.",
    )
    parser.add_argument(
        "--no-inchikey-lookup",
        action="store_true",
        help="Use only ChEMBL IDs already present in the workbook.",
    )
    parser.add_argument(
        "--no-target-lookup",
        action="store_true",
        help="Use target names/IDs from ChEMBL activity rows without extra target metadata requests.",
    )
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    if args.chembl_bioactivity:
        row_count = create_chembl_bioactivity_tsv(
            Path(args.excel),
            Path(args.chembl_output),
            args.sheet,
            max_compounds=args.max_compounds,
            max_activities=args.max_activities,
            fallback_to_inchikey_lookup=not args.no_inchikey_lookup,
            resolve_targets=not args.no_target_lookup,
        )
        print(f"Wrote {args.chembl_output}: {row_count} rows")
        return

    counts = convert_excel_to_tsv(Path(args.excel), Path(args.output_dir), args.sheet)

    for filename, row_count in counts.items():
        print(f"Wrote {filename}: {row_count} rows")


if __name__ == "__main__":
    main()
