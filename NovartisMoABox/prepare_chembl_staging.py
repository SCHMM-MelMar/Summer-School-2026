import csv
from pathlib import Path


ROOT = Path(__file__).resolve().parent.parent
STAGING = ROOT / "staging" / "custom_data"
TARGET_TSV = STAGING / "target.tsv"
BIOACTIVITY_TSV = STAGING / "bioactivity.tsv"
ADDED_TARGETS_TSV = STAGING / "target_added_from_bioactivity.tsv"


def read_rows(path):
    with path.open(newline="") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def write_rows(path, rows, fieldnames):
    with path.open("w", newline="") as handle:
        writer = csv.DictWriter(
            handle,
            fieldnames=fieldnames,
            delimiter="\t",
            lineterminator="\n",
        )
        writer.writeheader()
        writer.writerows(rows)


def main():
    targets = read_rows(TARGET_TSV)
    bioactivities = read_rows(BIOACTIVITY_TSV)

    known_targets = {row["target_key"] for row in targets}
    bioactivity_targets = {row["target_key"] for row in bioactivities}
    missing_targets = sorted(bioactivity_targets - known_targets)

    added_targets = [
        {"target_key": key, "type": "protein", "name": key} for key in missing_targets
    ]
    write_rows(TARGET_TSV, targets + added_targets, ["target_key", "type", "name"])
    write_rows(ADDED_TARGETS_TSV, added_targets, ["target_key", "type", "name"])

    for row in bioactivities:
        if not row.get("relation", "").strip():
            row["relation"] = "="
        if not row.get("unit", "").strip():
            row["unit"] = "unspecified"

    write_rows(
        BIOACTIVITY_TSV,
        bioactivities,
        [
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
        ],
    )

    print(f"added targets: {len(added_targets)}")
    print(f"normalized bioactivity rows: {len(bioactivities)}")
    print(f"added target report: {ADDED_TARGETS_TSV}")


if __name__ == "__main__":
    main()
