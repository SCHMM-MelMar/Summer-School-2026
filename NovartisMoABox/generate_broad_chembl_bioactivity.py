from __future__ import annotations

import argparse
import csv
from concurrent.futures import ThreadPoolExecutor, as_completed
from pathlib import Path

import pandas as pd

from NovartisDataset import (
    BIOACTIVITY_COLUMNS,
    DEFAULT_EXCEL,
    DEFAULT_SHEET,
    clean_value,
    chembl_get,
    first_present_dict,
    keep_chembl_activity,
    make_workbook_compound_maps,
)


ROOT = Path(__file__).resolve().parent
DATA_DIR = ROOT / "custom_data"
COMPOUND_TSV = DATA_DIR / "compound.tsv"
BIOACTIVITY_TSV = DATA_DIR / "bioactivity.tsv"
ASSAY_TYPES = "AC50,Activity,CC50,EC50,GI50,IC50,Inhibition,Kd,Ki,MIC,Potency"
RELATION_MAP = {">>": ">", "<<": "<"}


def read_compounds(path: Path) -> list[tuple[str, str]]:
    pairs = []
    seen = set()
    with path.open(newline="", encoding="utf-8") as handle:
        for row in csv.DictReader(handle, delimiter="\t"):
            inchikey = clean_value(row.get("inchikey"))
            chembl_ids = clean_value(row.get("chembl_id")).replace(";", "|").split("|")
            for chembl_id in chembl_ids:
                chembl_id = chembl_id.strip()
                key = (chembl_id, inchikey)
                if chembl_id and inchikey and key not in seen:
                    pairs.append(key)
                    seen.add(key)
    return sorted(pairs)


def load_moa_by_inchikey(excel_path: Path | None, sheet_name: str) -> dict[str, str]:
    if not excel_path or not excel_path.exists():
        return {}
    df = pd.read_excel(excel_path, sheet_name=sheet_name)
    _, inchikey_to_moa = make_workbook_compound_maps(df)
    return inchikey_to_moa


def records_for_compound(
    chembl_id: str,
    inchikey: str,
    moa: str,
    max_rows: int,
) -> list[dict[str, str]]:
    records = []
    offset = 0
    limit = 100

    while len(records) < max_rows:
        payload = chembl_get(
            "activity",
            {
                "molecule_chembl_id": chembl_id,
                "standard_type__in": ASSAY_TYPES,
                "limit": limit,
                "offset": offset,
            },
        )
        activities = payload.get("activities", [])
        if not activities:
            break

        for activity in activities:
            if not keep_chembl_activity(activity):
                continue
            records.append(
                {
                    "inchikey": inchikey,
                    "target_key": first_present_dict(
                        activity, ["target_pref_name", "target_chembl_id"]
                    ),
                    "moa": moa,
                    "bioactivity_type": clean_value(activity.get("standard_type")),
                    "relation": RELATION_MAP.get(
                        clean_value(activity.get("standard_relation")),
                        clean_value(activity.get("standard_relation")) or "=",
                    ),
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
            if len(records) >= max_rows:
                break

        page_meta = payload.get("page_meta", {})
        if not page_meta.get("next") or len(activities) < limit:
            break
        offset += limit

    return records


def write_rows(path: Path, rows: list[dict[str, str]]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(
            handle,
            fieldnames=BIOACTIVITY_COLUMNS,
            delimiter="\t",
            lineterminator="\n",
        )
        writer.writeheader()
        writer.writerows(rows)


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Build broad ChEMBL bioactivity.tsv coverage for custom_data/compound.tsv."
    )
    parser.add_argument("--compound-tsv", default=str(COMPOUND_TSV))
    parser.add_argument("--output", default=str(BIOACTIVITY_TSV))
    parser.add_argument("--excel", default=DEFAULT_EXCEL)
    parser.add_argument("--sheet", default=DEFAULT_SHEET)
    parser.add_argument("--max-rows-per-chembl-id", type=int, default=5)
    parser.add_argument("--workers", type=int, default=12)
    parser.add_argument("--max-chembl-ids", type=int, default=None)
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    compounds = read_compounds(Path(args.compound_tsv))
    if args.max_chembl_ids is not None:
        compounds = compounds[: args.max_chembl_ids]

    moa_by_inchikey = load_moa_by_inchikey(Path(args.excel), args.sheet)
    rows = []
    completed = 0
    with ThreadPoolExecutor(max_workers=args.workers) as pool:
        futures = {
            pool.submit(
                records_for_compound,
                chembl_id,
                inchikey,
                moa_by_inchikey.get(inchikey, ""),
                args.max_rows_per_chembl_id,
            ): (chembl_id, inchikey)
            for chembl_id, inchikey in compounds
        }
        for future in as_completed(futures):
            completed += 1
            chembl_id, inchikey = futures[future]
            try:
                rows.extend(future.result())
            except Exception as error:
                print(f"Skipping {chembl_id}/{inchikey}: {error}", flush=True)
            if completed % 250 == 0 or completed == len(compounds):
                covered = len({row["inchikey"] for row in rows})
                print(
                    f"Queried {completed}/{len(compounds)} ChEMBL IDs; "
                    f"kept {len(rows)} rows for {covered} InChIKeys...",
                    flush=True,
                )

    rows.sort(
        key=lambda row: (
            row["inchikey"],
            row["target_key"],
            row["bioactivity_type"],
            row["source"],
            row["source_xref"],
        )
    )
    write_rows(Path(args.output), rows)
    print(f"Wrote {args.output}: {len(rows)} rows")
    print(f"Unique InChIKeys with rows: {len({row['inchikey'] for row in rows})}")


if __name__ == "__main__":
    main()
