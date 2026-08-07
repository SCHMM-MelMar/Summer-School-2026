"""Structural similarity between compounds.

Nothing is precomputed and stored. Fingerprinting all 16564 structures takes
3.4 seconds and the largest target in the database, 3115 compounds, is 4.85
million pairs in 0.65 seconds. A stored similarity matrix would be 137 million
pairs, larger than the whole database, to answer questions almost none of which
get asked. The fingerprints are cached per database so the 3.4 seconds happens
once per session.

RDKit is an optional dependency. Everything else in probedb works without it.
"""

import pandas as pd

from .db import ProbeDB

RADIUS = 2
BITS = 2048

# Count fingerprints, not bits. A binary fingerprint records which atom
# environments a molecule has and not how many, so lauric acid (C12) and
# behenic acid (C22) are identical to it, as are azelaic and sebacic acid and a
# piperidine and its pyrrolidine analogue. Counts separate them: 0.57, 0.92,
# 0.95. It costs speed, measured on 4.5 million pairs: 0.42s on bits, 5.8s on
# counts. Set COUNTS = False to trade the accuracy back for the 14x.
COUNTS = True

_CACHE = {}


def _rdkit():
    try:
        from rdkit import Chem, DataStructs, RDLogger
        from rdkit.Chem import rdFingerprintGenerator
    except ImportError as missing:  # pragma: no cover
        raise ImportError(
            "structural similarity needs RDKit: `pip install rdkit`, or "
            "`mamba install -c conda-forge rdkit`"
        ) from missing
    RDLogger.DisableLog("rdApp.*")
    return Chem, DataStructs, rdFingerprintGenerator


def fingerprints(db):
    """Morgan fingerprints for every compound with a structure.

    Chirality is on. Without it a probe and its inactive enantiomer score 1.00,
    and that is the one pair in this database that most needs telling apart:
    (+)-JQ1 and (-)-JQ1 come to 0.89 with it and 1.00 without.

    Cached per connection, so the 3.4 seconds it takes for all 16564 structures
    happens once a session and not once a question.
    """
    key = id(db.conn)
    if key in _CACHE:
        return _CACHE[key]

    Chem, _, rdFingerprintGenerator = _rdkit()
    generator = rdFingerprintGenerator.GetMorganGenerator(
        radius=RADIUS, fpSize=BITS, includeChirality=True
    )
    fingerprint = generator.GetCountFingerprint if COUNTS else generator.GetFingerprint

    rows = db.read("SELECT inchikey, smiles FROM compound WHERE smiles IS NOT NULL")
    out = {}
    for inchikey, smiles in zip(rows.inchikey, rows.smiles):
        mol = Chem.MolFromSmiles(smiles)
        if mol is not None:
            out[inchikey] = fingerprint(mol)
    _CACHE[key] = out
    return out


def _keys_for(db, target=None, family=None, set=None, compounds=None):
    if compounds is not None:
        return [db.compound_key(c) for c in compounds]
    if target is not None:
        ids = db.targets_for(target)
        marks = ",".join("?" * len(ids))
        rows = db.conn.execute(
            f"SELECT DISTINCT inchikey FROM bioactivity WHERE target_id IN ({marks})",
            ids,
        ).fetchall()
        return [r[0] for r in rows]
    if family is not None:
        return db.compounds(family=family)["inchikey"].tolist()
    if set is not None:
        return db.compounds(set=set)["inchikey"].tolist()
    return db.read("SELECT inchikey FROM compound")["inchikey"].tolist()


def _selected(db, target, family, set, compounds):
    fps = fingerprints(db)
    keys = dict.fromkeys(_keys_for(db, target, family, set, compounds))
    return [k for k in keys if k in fps], fps


def _rows(db, keys):
    """Every pair as (i, j, score), one row per pair, nothing materialised."""
    import numpy as np

    _, DataStructs, _ = _rdkit()
    fps = fingerprints(db)
    vectors = [fps[k] for k in keys]
    for i in range(len(keys) - 1):
        yield i, np.asarray(
            DataStructs.BulkTanimotoSimilarity(vectors[i], vectors[i + 1 :])
        )


def _names(db, keys):
    if not keys:
        return {}
    return dict(
        db.conn.execute(
            "SELECT inchikey, name FROM compound WHERE inchikey IN ({})".format(
                ",".join("?" * len(keys))
            ),
            keys,
        )
    )


def pairs(db, target=None, family=None, set=None, compounds=None, threshold=0.0):
    """Every pair of compounds and how similar they are, most similar first.

    Pick the set with one of target=, family=, set= or compounds=. A pair is
    one row, not two, and a compound is never compared with itself.

    Note the threshold. A target with 3000 compounds is five million pairs:
    the arithmetic takes under a second, building a DataFrame row for each
    takes twenty. Ask for what you are going to read.
    """
    keys, _ = _selected(db, target, family, set, compounds)
    names = _names(db, keys)

    rows = []
    for i, scores in _rows(db, keys):
        for offset in (scores >= threshold).nonzero()[0]:
            j = i + 1 + int(offset)
            rows.append(
                (keys[i], names.get(keys[i]), keys[j], names.get(keys[j]),
                 float(scores[offset]),
                 # the first block of an InChIKey hashes the skeleton, so equal
                 # first blocks and different keys means the two differ only in
                 # stereochemistry or isotopes. a fingerprint cannot always see
                 # that difference, and this says so exactly rather than
                 # leaving a 1.00 looking like a duplicate
                 keys[i][:14] == keys[j][:14])
            )

    out = pd.DataFrame(
        rows,
        columns=["inchikey_a", "compound_a", "inchikey_b", "compound_b",
                 "tanimoto", "same_skeleton"],
    )
    return out.sort_values("tanimoto", ascending=False).reset_index(drop=True)


MATRIX_LIMIT = 400


def matrix(db, target=None, family=None, set=None, compounds=None, limit=MATRIX_LIMIT):
    """The same thing square, for a set small enough to read as a heatmap."""
    _, DataStructs, _ = _rdkit()
    fps = fingerprints(db)

    keys, _ = _selected(db, target, family, set, compounds)
    if len(keys) > limit:
        raise ValueError(
            f"{len(keys)} compounds is {len(keys) ** 2} cells, which is not a "
            f"heatmap anybody can read. Narrow the selection, raise limit= if "
            f"you mean it, or use pairs() with a threshold instead"
        )
    names = _names(db, keys)
    labels = [names.get(k) or k for k in keys]
    vectors = [fps[k] for k in keys]
    data = [DataStructs.BulkTanimotoSimilarity(v, vectors) for v in vectors]
    return pd.DataFrame(data, index=labels, columns=labels)


def neighbours(db, compound, n=10, threshold=0.0):
    """The compounds most like this one, across everything with a structure."""
    _, DataStructs, _ = _rdkit()
    fps = fingerprints(db)

    key = db.compound_key(compound)
    if key not in fps:
        raise KeyError(f"no structure on file for {compound!r}")

    keys = [k for k in fps if k != key]
    scores = DataStructs.BulkTanimotoSimilarity(fps[key], [fps[k] for k in keys])
    out = pd.DataFrame({"inchikey": keys, "tanimoto": scores})
    out = out[out.tanimoto >= threshold].nlargest(n, "tanimoto")
    named = db.read(
        "SELECT inchikey, name AS compound FROM compound WHERE inchikey IN ({})".format(
            ",".join("?" * len(out))
        ),
        *out.inchikey,
    )
    return out.merge(named, on="inchikey").sort_values(
        "tanimoto", ascending=False
    ).reset_index(drop=True)


def summary(db, target=None, family=None, set=None, compounds=None, threshold=0.7):
    """One row: is this a diverse set of compounds or a series of analogues?

    Reads the scores as numbers and never builds a row per pair, so it stays
    under a second on the largest target in the database.
    """
    import numpy as np

    keys, _ = _selected(db, target, family, set, compounds)
    if len(keys) < 2:
        return pd.DataFrame()

    # a dict standing in for a set, because `set` is a parameter name here
    chunks, close = [], {}
    for i, scores in _rows(db, keys):
        chunks.append(scores)
        for offset in (scores >= threshold).nonzero()[0]:
            close[keys[i]] = close[keys[i + 1 + int(offset)]] = True

    everything = np.concatenate(chunks)
    return pd.DataFrame(
        [
            {
                "compounds": len(keys),
                "pairs": everything.size,
                "median": float(np.median(everything)),
                "max": float(everything.max()),
                f"pairs_over_{threshold}": int((everything >= threshold).sum()),
                # a compound with a close relative in the set is part of a
                # series; one with none was made on its own
                "in_a_series": len(close),
            }
        ]
    )
