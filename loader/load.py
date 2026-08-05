from pathlib import Path

from .validate import read, validate

STAGING = Path(__file__).resolve().parent.parent / "staging"


def load(db, directory, source=None, strict=True):
    directory = Path(directory)
    problems = validate(directory)
    errors = [p for p in problems if not p.endswith(("ignored", "kept anyway"))]
    if errors and strict:
        raise ValueError(f"{directory} did not validate:\n  " + "\n  ".join(errors))

    report = {
        "directory": str(directory),
        "problems": problems,
        "compounds": 0,
        "targets": 0,
        "complexes": 0,
        "bioactivities": 0,
        "duplicates_skipped": 0,
        "chembl_ids": 0,
    }
    members = {}
    for row in read(directory, "uniprot"):
        members.setdefault(row["target_key"], []).append(
            (
                row["uniprot_id"],
                row.get("hgnc"),
                row.get("species"),
                row.get("entrez_gene"),
            )
        )

    targets = {}
    for row in read(directory, "target"):
        key = row["target_key"]
        targets[key] = db.add_target(
            row.get("type") or "protein", row.get("name") or key, members.get(key, [])
        )
        report["targets"] += 1
        report["complexes"] += len(members.get(key, [])) > 1

    compounds = {}
    for row in read(directory, "compound"):
        inchikey = row["inchikey"].upper()
        compounds[inchikey] = db.add_compound(
            inchikey, row.get("smiles"), row.get("name")
        )
        report["compounds"] += 1
        for chembl_id in (row.get("chembl_id") or "").replace(";", "|").split("|"):
            if chembl_id.strip():
                db.add_chembl(inchikey, chembl_id)
                report["chembl_ids"] += 1

    seen = {
        tuple(r)
        for r in db.conn.execute(
            "SELECT inchikey, target_id, source_id, bioactivity_type, relation, "
            "       value, unit FROM bioactivity"
        )
    }

    sources = {}
    for row in read(directory, "bioactivity"):
        source_key = (row.get("source_db") or source, row.get("source") or None)
        if source_key not in sources:
            sources[source_key] = db.add_source(*source_key)

        key = (
            compounds[row["inchikey"].upper()],
            targets[row["target_key"]],
            sources[source_key],
            row.get("bioactivity_type") or None,
            row.get("relation") or None,
            None if not row.get("value") else float(row["value"]),
            row.get("unit") or None,
        )
        if key in seen:
            report["duplicates_skipped"] += 1
            continue
        seen.add(key)

        db.add_bioactivity(
            compounds[row["inchikey"].upper()],
            targets[row["target_key"]],
            moa=row.get("moa"),
            bioactivity_type=row.get("bioactivity_type") or None,
            relation=row.get("relation"),
            value=row["value"] or None,
            unit=row["unit"],
            assay_type=row.get("assay_type"),
            assay_description=row.get("assay_description"),
            cell_line=row.get("cell_line"),
            concentration=row.get("concentration"),
            concentration_unit=row.get("concentration_unit"),
            source_id=sources[source_key],
            source_xref=row.get("source_xref"),
        )
        report["bioactivities"] += 1

    db.commit()
    return report


def load_all(db, root=STAGING, strict=True):
    return [
        load(db, path, source=path.name, strict=strict)
        for path in sorted(Path(root).iterdir())
        if path.is_dir() and not path.name.startswith("_")
    ]
