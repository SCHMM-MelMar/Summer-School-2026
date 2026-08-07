from pathlib import Path

from probedb.db import INCHIKEY

from .validate import TARGET_TYPES, UNIPROT, read, validate

STAGING = Path(__file__).resolve().parent.parent / "staging"


def load(db, directory, source=None, strict=True):
    directory = Path(directory)
    problems = validate(directory)
    errors = [
        p for p in problems
        if not p.endswith(("ignored", "kept anyway", "skipped"))
    ]
    if errors and strict:
        # a directory can be wrong thousands of times over and printing all of
        # it buries the one line that says what to fix
        shown = "\n  ".join(errors[:10])
        rest = len(errors) - 10
        raise ValueError(
            f"{directory} did not validate, {len(errors)} problem(s):\n  {shown}"
            + (f"\n  ... and {rest} more, run validate() to see them all" if rest > 0 else "")
        )

    report = {
        "directory": str(directory),
        "problems": problems,
        "compounds": 0,
        "targets": 0,
        "complexes": 0,
        "sets": 0,
        "set_members": 0,
        "bioactivities": 0,
        "duplicates_skipped": 0,
        "chembl_ids": 0,
        "targets_skipped": 0,
        "compounds_skipped": 0,
        "accessions_skipped": 0,
        "bioactivities_skipped": 0,
    }
    members = {}
    for row in read(directory, "uniprot"):
        # a target whose accession is unusable keeps its row and is identified
        # by its name instead, which is what add_target falls back to
        if not UNIPROT.match(row["uniprot_id"].upper()):
            report["accessions_skipped"] += 1
            continue
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
        kind = row.get("type") or "protein"
        # strict mode already refused this; under strict=False the row is
        # dropped rather than taking the whole load down with it
        if kind not in TARGET_TYPES:
            report["targets_skipped"] += 1
            continue
        targets[key] = db.add_target(
            kind, row.get("name") or key, members.get(key, [])
        )
        report["targets"] += 1
        report["complexes"] += len(members.get(key, [])) > 1

    compounds = {}
    for row in read(directory, "compound"):
        inchikey = row["inchikey"].upper()
        if not INCHIKEY.match(inchikey):
            report["compounds_skipped"] += 1
            continue
        compounds[inchikey] = db.add_compound(
            inchikey, row.get("smiles"), row.get("name")
        )
        report["compounds"] += 1
        for chembl_id in (row.get("chembl_id") or "").replace(";", "|").split("|"):
            if chembl_id.strip():
                db.add_chembl(inchikey, chembl_id)
                report["chembl_ids"] += 1

    # one directory is one set, named after itself. a compound arriving from
    # three directories is one compound row and three memberships, which is
    # how "where did this come from" stays answerable after the merge
    if compounds:
        set_id = db.add_set(source or directory.name, "library", source)
        report["sets"] = 1
        for inchikey in compounds.values():
            db.add_set_member(set_id, inchikey)
            report["set_members"] += 1

    seen = {
        tuple(r)
        for r in db.conn.execute(
            "SELECT inchikey, target_id, source_id, bioactivity_type, relation, "
            "       value, unit FROM bioactivity"
        )
    }

    sources, resolved = {}, set()
    for row in read(directory, "bioactivity"):
        source_key = (row.get("source_db") or source, row.get("source") or None)
        xref_id = row.get("xref_id")
        # a source may only carry its xref prefix on some of its rows
        if source_key not in sources or (xref_id and source_key not in resolved):
            sources[source_key] = db.add_source(*source_key, xref_id)
            if xref_id:
                resolved.add(source_key)

        # a measurement whose compound or target never made it in has nothing
        # to attach to
        inchikey = compounds.get(row["inchikey"].upper())
        target_id = targets.get(row["target_key"])
        if inchikey is None or target_id is None:
            report["bioactivities_skipped"] += 1
            continue

        key = (
            inchikey,
            target_id,
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
            inchikey,
            target_id,
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
