import csv
import re
import sys
import urllib.parse
import urllib.request
from collections import defaultdict
from pathlib import Path


ROOT = Path(__file__).resolve().parent.parent
TARGET_TSV = ROOT / "staging" / "custom_data" / "target.tsv"
UNIPROT_TSV = ROOT / "staging" / "custom_data" / "uniprot.tsv"
UNMATCHED_TSV = ROOT / "staging" / "custom_data" / "uniprot_unmatched.tsv"
UNIPROT_CACHE = ROOT / "staging" / "custom_data" / "_uniprot_reviewed_human.tsv"

MANUAL_MAPPINGS = {
    "CKBE": ("P12277", "CKB", "Homo sapiens (Human)"),
    "MMP23A": ("O75900", "MMP23B", "Homo sapiens (Human)"),
    "pbp2b": ("P0A3M6", "PBP2B", "Streptococcus pneumoniae (strain ATCC BAA-255 / R6)"),
    "rpoB": ("P9WGY9", "RPOB", "Mycobacterium tuberculosis (strain ATCC 25618 / H37Rv)"),
}


def read_targets():
    with TARGET_TSV.open(newline="") as handle:
        rows = list(csv.DictReader(handle, delimiter="\t"))
    return [row for row in rows if (row.get("type") or "protein") == "protein"]


def normalize_symbol(value):
    return re.sub(r"[^A-Z0-9]", "", (value or "").upper())


def fetch_reviewed_human():
    if UNIPROT_CACHE.exists():
        with UNIPROT_CACHE.open(newline="") as handle:
            yield from csv.DictReader(handle, delimiter="\t")
        return

    query = urllib.parse.urlencode(
        {
            "query": "reviewed:true AND organism_id:9606",
            "fields": "accession,gene_primary,gene_names,organism_name",
            "format": "tsv",
            "size": "500",
        }
    )
    url = f"https://rest.uniprot.org/uniprotkb/search?{query}"

    lines = []
    while url:
        request = urllib.request.Request(url, headers={"User-Agent": "summer-school-db/0.1"})
        with urllib.request.urlopen(request, timeout=60) as response:
            text = response.read().decode("utf-8")
            page_lines = text.splitlines()
            lines.extend(page_lines if not lines else page_lines[1:])
            yield from csv.DictReader(text.splitlines(), delimiter="\t")
            link = response.headers.get("Link", "")
        match = re.search(r"<([^>]+)>; rel=\"next\"", link)
        url = match.group(1) if match else None

    UNIPROT_CACHE.write_text("\n".join(lines) + "\n")


def build_index(entries):
    index = defaultdict(list)
    for entry in entries:
        accession = entry["Entry"].strip()
        primary = entry["Gene Names (primary)"].strip()
        names = entry["Gene Names"].replace(";", " ").split()
        organism = entry["Organism"].strip()
        if not accession or not primary:
            continue
        for symbol in {primary, *names}:
            key = normalize_symbol(symbol)
            if key:
                index[key].append((accession, primary, organism))
    return index


def choose_match(target_key, matches):
    if not matches:
        return None
    exact = [m for m in matches if m[1].upper() == target_key.upper()]
    candidates = exact or matches
    return sorted(candidates, key=lambda m: (m[1].upper() != target_key.upper(), m[0]))[0]


def main():
    targets = read_targets()
    index = build_index(fetch_reviewed_human())

    matched, unmatched = [], []
    seen = set()
    for target in targets:
        target_key = target["target_key"].strip()
        match = MANUAL_MAPPINGS.get(
            target_key, choose_match(target_key, index.get(normalize_symbol(target_key), []))
        )
        if match:
            accession, hgnc, species = match
            row_key = (accession, target_key)
            if row_key not in seen:
                matched.append(
                    {
                        "uniprot_id": accession,
                        "target_key": target_key,
                        "hgnc": hgnc,
                        "species": species,
                    }
                )
                seen.add(row_key)
        else:
            unmatched.append(target)

    with UNIPROT_TSV.open("w", newline="") as handle:
        writer = csv.DictWriter(
            handle,
            fieldnames=["uniprot_id", "target_key", "hgnc", "species"],
            delimiter="\t",
            lineterminator="\n",
        )
        writer.writeheader()
        writer.writerows(matched)

    with UNMATCHED_TSV.open("w", newline="") as handle:
        writer = csv.DictWriter(
            handle,
            fieldnames=["target_key", "type", "name"],
            delimiter="\t",
            lineterminator="\n",
        )
        writer.writeheader()
        writer.writerows(unmatched)

    print(f"matched protein targets: {len(matched)}")
    print(f"unmatched protein targets: {len(unmatched)}")
    print(f"wrote: {UNIPROT_TSV}")
    print(f"unmatched report: {UNMATCHED_TSV}")
    return 0 if matched else 1


if __name__ == "__main__":
    raise SystemExit(main())
