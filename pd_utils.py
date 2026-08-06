"""
pd_utils.py — Reusable query helpers for the Probes & Drugs (P&D) database.

Designed for computational biologists: each function wraps a common query
into a one-liner that returns a pandas DataFrame. Raw SQL is visible inside
each function for those who want to learn or modify.

Usage:
    from pd_utils import *
    df = get_compound_targets("JQ1")
    df = get_target_compounds("EGFR")
    df = get_family_compounds(["KRAS", "NRAS", "HRAS"])

Database: 25-table simplified P&D SQLite database.
"""

import sqlite3
import pandas as pd
import matplotlib
import matplotlib.pyplot as plt
import numpy as np
from collections import Counter

# ── Configuration ────────────────────────────────────────────────
DB_PATH = '/mnt/results/pd_simple.sqlite'

# Colorblind-friendly palette (Phylo)
PALETTE = ['#0279EE', '#FF9400', '#75A025', '#FD9BED', '#E9ED4C', '#666666']

# Matplotlib defaults
matplotlib.rcParams['font.family'] = ['Liberation Sans', 'Arimo', 'DejaVu Sans']
matplotlib.rcParams['svg.fonttype'] = 'none'
matplotlib.rcParams['figure.dpi'] = 120

# Log-scale activity types (use these for potency ranking — not percentage inhibition)
LOG_ACTIVITY_TYPES = ('pIC50', 'pKd', 'pKi', 'pEC50', 'pAC50', 'pPotency')

# Chemogenomic set keywords (no dedicated compoundsettype label)
CHEMOGEN_KEYWORDS = ['chemogen', 'MoA', 'LINCS', 'JUMP', 'PKIS', 'Informer', 'LSP', 'NIBR']


# ── Database connection ──────────────────────────────────────────
_conn = None

def get_connection():
    """Return a shared read-only SQLite connection."""
    global _conn
    if _conn is None:
        _conn = sqlite3.connect(f'file:{DB_PATH}?mode=ro', uri=True)
    return _conn


def run_query(sql, params=None, show_raw=True, max_raw=15):
    """
    Execute SQL, optionally print raw cursor output, return pandas DataFrame.

    Parameters
    ----------
    sql : str           — SQL query string
    params : tuple/list — bind parameters for ? placeholders
    show_raw : bool     — if True, print raw tabular output
    max_raw : int       — max rows to print in raw output

    Returns
    -------
    pd.DataFrame
    """
    cur = get_connection().execute(sql, params or ())
    cols = [d[0] for d in cur.description]
    rows = cur.fetchall()

    if show_raw and rows:
        print(f'Raw SQL output ({len(rows)} rows):')
        print('  ' + ' | '.join(cols))
        print('  ' + '-' * 80)
        for r in rows[:max_raw]:
            print('  ' + ' | '.join(str(v)[:25] if v is not None else 'NULL' for v in r))
        if len(rows) > max_raw:
            print(f'  ... ({len(rows) - max_raw} more rows)')
        print()
    elif show_raw:
        print('Raw SQL output: 0 rows\n')

    return pd.DataFrame(rows, columns=cols)


def list_tables():
    """List all tables and their row counts in the database."""
    conn = get_connection()
    tables = [r[0] for r in conn.execute(
        "SELECT name FROM sqlite_master WHERE type='table' ORDER BY name")]
    rows = [(t, conn.execute(f'SELECT COUNT(*) FROM "{t}"').fetchone()[0]) for t in tables]
    return pd.DataFrame(rows, columns=['table', 'n_rows'])


# ════════════════════════════════════════════════════════════════
#  COMPOUND-CENTRIC QUERIES (US1)
# ════════════════════════════════════════════════════════════════

def get_compound(name):
    """
    Look up a compound by name (case-insensitive).

    Returns DataFrame with compoundid, pdid, name, smiles, inchikey.
    """
    sql = """
    SELECT compoundid, pdid, name, smiles, inchikey
    FROM compound
    WHERE name LIKE ? COLLATE NOCASE
    """
    return run_query(sql, params=(f'%{name}%',), show_raw=False)


def get_compound_targets(name, inactive=False, log_only=True,
                         exact_only=True, min_confidence=1,
                         flag_contradictions=True):
    """
    Get all targets a compound has activity against, with best potency per target.

    Path: compound → activity → target → targettobasetarget → basetarget

    Parameters
    ----------
    name : str            — compound name (case-insensitive)
    inactive : bool       — if True, include inactive measurements
    log_only : bool       — if True (default), restrict to log-scale activity types
                            (pIC50, pKd, pKi, pEC50, pAC50, pPotency) to avoid mixing
                            with percentage-scale values (Inhibition, Dmax)
    exact_only : bool     — if True (default), exclude targets whose ONLY measurements
                            use the '>' operator (screening negatives like "IC50 > 30 µM").
                            These are not real target engagements — they mean the compound
                            was tested and found to be inactive above a threshold.
    min_confidence : int  — minimum confidence level (1=high/direct, 2=lower/derived).
                            Default 1 filters to directly measured values only.
    flag_contradictions : bool — if True (default), add a 'contradiction_flag' column
                            that marks targets where binding (pKd) and functional (pIC50)
                            measurements disagree by >2 log units, indicating possible
                            data errors.

    Returns DataFrame with pdid, compound_name, target_gene, target_name,
             best_potency, best_activity_type, activity_types, n_measurements,
             contradiction_flag (if flag_contradictions=True).
    """
    inactive_filter = '' if inactive else 'AND a.inactive = 0'
    if log_only:
        quoted = ",".join(f'"{t}"' for t in LOG_ACTIVITY_TYPES)
        type_filter = f'AND a.activity_type IN ({quoted})'
    else:
        type_filter = ''
    conf_filter = f'AND a.confidence <= {int(min_confidence)}'

    # Look up compound_id first (enables index usage on activity.compound_id)
    conn = get_connection()
    cid_row = conn.execute(
        "SELECT compoundid FROM compound WHERE name LIKE ? COLLATE NOCASE",
        (name,)
    ).fetchone()

    if cid_row is None:
        return pd.DataFrame()

    compound_id = cid_row[0]

    sql = f"""
    SELECT
        c.pdid,
        c.name AS compound_name,
        bt.gene_name AS target_gene,
        bt.name AS target_name,
        MAX(a.activity_value) AS best_potency,
        GROUP_CONCAT(DISTINCT a.activity_type) AS activity_types,
        GROUP_CONCAT(DISTINCT a.value_type) AS value_types,
        COUNT(*) AS n_measurements
    FROM activity a
    JOIN compound c ON a.compound_id = c.compoundid
    JOIN target t ON a.target_id = t.targetid
    JOIN targettobasetarget ttb ON t.targetid = ttb.target_id
    JOIN basetarget bt ON ttb.basetarget_id = bt.basetargetid
    WHERE a.compound_id = ?
      {inactive_filter}
      {type_filter}
      {conf_filter}
      AND a.activity_value IS NOT NULL
    GROUP BY c.pdid, c.name, bt.gene_name, bt.name
    ORDER BY best_potency DESC
    """
    df = run_query(sql, params=(compound_id,))

    # Filter out targets whose ONLY measurements are '>' (screening negatives)
    if exact_only and not df.empty:
        has_exact = df['value_types'].str.contains('=', na=False)
        df = df[has_exact].copy()

    # Add the activity type of the best measurement using a single batch query
    if not df.empty:
        df['best_activity_type'] = _batch_best_activity_types(df, name,
                                                              log_only, min_confidence)

    # Flag contradictions (binding vs functional disagree by >2 log units)
    if flag_contradictions and not df.empty:
        df['contradiction_flag'] = _batch_flag_contradictions(df, name)

    return df


def get_compound_targets_detailed(name, inactive=False, log_only=True,
                                  exact_only=True, min_confidence=1):
    """
    Get per-target per-activity-type breakdown for a compound.

    Returns one row per (target, activity_type) with the best value for that
    combination. Use this to plot potency grouped/hued by measurement type.

    Parameters
    ----------
    name : str            — compound name (case-insensitive)
    inactive : bool       — if True, include inactive measurements
    log_only : bool       — if True (default), restrict to log-scale activity types
    exact_only : bool     — if True (default), exclude targets whose ONLY measurements
                            use the '>' operator (screening negatives)
    min_confidence : int  — minimum confidence level (1=high, 2=lower). Default 1.

    Returns DataFrame with pdid, compound_name, target_gene, target_name,
             activity_type, best_potency, n_measurements.
    """
    inactive_filter = '' if inactive else 'AND a.inactive = 0'
    if log_only:
        quoted = ",".join(f'"{t}"' for t in LOG_ACTIVITY_TYPES)
        type_filter = f'AND a.activity_type IN ({quoted})'
    else:
        type_filter = ''
    conf_filter = f'AND a.confidence <= {int(min_confidence)}'

    # Look up compound_id first (enables index usage)
    conn = get_connection()
    cid_row = conn.execute(
        "SELECT compoundid FROM compound WHERE name LIKE ? COLLATE NOCASE",
        (name,)
    ).fetchone()
    if cid_row is None:
        return pd.DataFrame()
    compound_id = cid_row[0]

    sql = f"""
    SELECT
        c.pdid,
        c.name AS compound_name,
        bt.gene_name AS target_gene,
        bt.name AS target_name,
        a.activity_type,
        MAX(a.activity_value) AS best_potency,
        GROUP_CONCAT(DISTINCT a.value_type) AS value_types,
        COUNT(*) AS n_measurements
    FROM activity a
    JOIN compound c ON a.compound_id = c.compoundid
    JOIN target t ON a.target_id = t.targetid
    JOIN targettobasetarget ttb ON t.targetid = ttb.target_id
    JOIN basetarget bt ON ttb.basetarget_id = bt.basetargetid
    WHERE a.compound_id = ?
      {inactive_filter}
      {type_filter}
      {conf_filter}
      AND a.activity_value IS NOT NULL
    GROUP BY c.pdid, c.name, bt.gene_name, bt.name, a.activity_type
    ORDER BY bt.gene_name, best_potency DESC
    """
    df = run_query(sql, params=(compound_id,), show_raw=False)

    # Filter out targets whose ONLY measurements are '>' (screening negatives)
    if exact_only and not df.empty:
        # Keep only targets that have at least one exact '=' measurement
        # (identified at the target level, not per-type)
        target_has_exact = (
            df.groupby('target_gene')['value_types']
            .apply(lambda s: any('=' in str(v) for v in s))
        )
        valid_targets = target_has_exact[target_has_exact].index
        df = df[df['target_gene'].isin(valid_targets)].copy()

    return df


def _batch_best_activity_types(df, compound_name, log_only=True, min_confidence=1):
    """
    Batch helper: find which activity_type produced the best_potency for each
    target in df, using a single query instead of N per-row queries.
    """
    conn = get_connection()
    # Look up compound_id for index usage
    cid_row = conn.execute(
        "SELECT compoundid FROM compound WHERE name LIKE ? COLLATE NOCASE",
        (compound_name,)
    ).fetchone()
    if cid_row is None:
        return pd.Series(['?'] * len(df), index=df.index)
    compound_id = cid_row[0]

    if log_only:
        quoted = ",".join(f'"{t}"' for t in LOG_ACTIVITY_TYPES)
        type_filter = f'AND a.activity_type IN ({quoted})'
    else:
        type_filter = ''
    conf_filter = f'AND a.confidence <= {int(min_confidence)}'

    rows = conn.execute(f"""
        SELECT bt.gene_name, a.activity_value, a.activity_type
        FROM activity a
        JOIN target t ON a.target_id = t.targetid
        JOIN targettobasetarget ttb ON t.targetid = ttb.target_id
        JOIN basetarget bt ON ttb.basetarget_id = bt.basetargetid
        WHERE a.compound_id = ?
          AND a.inactive = 0
          AND a.activity_value IS NOT NULL
          {type_filter}
          {conf_filter}
    """, (compound_id,)).fetchall()

    lookup = {}
    for gene, val, atype in rows:
        key = (gene, val)
        if key not in lookup:
            lookup[key] = atype

    return df.apply(
        lambda row: lookup.get((row['target_gene'], row['best_potency']), '?'),
        axis=1
    )


def _batch_flag_contradictions(df, compound_name):
    """
    Batch helper: flag targets where binding (pKd) and functional (pIC50)
    measurements disagree by more than 2 log units, using a single query.
    """
    conn = get_connection()
    cid_row = conn.execute(
        "SELECT compoundid FROM compound WHERE name LIKE ? COLLATE NOCASE",
        (compound_name,)
    ).fetchone()
    if cid_row is None:
        return pd.Series([False] * len(df), index=df.index)
    compound_id = cid_row[0]

    rows = conn.execute("""
        SELECT
            bt.gene_name,
            MAX(CASE WHEN a.activity_type = 'pKd' AND a.value_type = '='
                      THEN a.activity_value END) AS best_pkd_exact,
            MAX(CASE WHEN a.activity_type = 'pIC50' AND a.value_type = '='
                      THEN a.activity_value END) AS best_pic50_exact,
            MAX(CASE WHEN a.activity_type = 'pKd' AND a.value_type = '>'
                      THEN a.activity_value END) AS best_pkd_gt
        FROM activity a
        JOIN target t ON a.target_id = t.targetid
        JOIN targettobasetarget ttb ON t.targetid = ttb.target_id
        JOIN basetarget bt ON ttb.basetarget_id = bt.basetargetid
        WHERE a.compound_id = ?
          AND a.inactive = 0
          AND a.activity_value IS NOT NULL
          AND a.activity_type IN ('pKd', 'pIC50')
        GROUP BY bt.gene_name
    """, (compound_id,)).fetchall()

    lookup = {}
    for gene, pkd_ex, pic50_ex, pkd_gt in rows:
        lookup[gene] = (pkd_ex, pic50_ex, pkd_gt)

    flags = []
    for _, row in df.iterrows():
        gene = row['target_gene']
        pkd_ex, pic50_ex, pkd_gt = lookup.get(gene, (None, None, None))

        is_contradictory = False
        if pkd_ex is not None and pic50_ex is not None:
            if abs(pkd_ex - pic50_ex) > 2.0:
                is_contradictory = True
        if pkd_gt is not None and pkd_ex is None and pic50_ex is not None:
            if (pic50_ex - pkd_gt) > 2.0:
                is_contradictory = True

        flags.append(is_contradictory)
    return pd.Series(flags, index=df.index)


def get_compound_actions(name):
    """
    Get curated mechanism-of-action annotations for a compound.

    Path: compound → compoundaction → target → targettobasetarget → basetarget

    Returns DataFrame with pdid, compound_name, target_gene, action,
             primary_target, drug_target.
    """
    sql = """
    SELECT
        c.pdid,
        c.name AS compound_name,
        bt.gene_name AS target_gene,
        ca.actiontype_id AS action,
        at.type AS action_type,
        ca.primary_target,
        ca.drug_target
    FROM compound c
    JOIN compoundaction ca ON c.compoundid = ca.compound_id
    JOIN target t ON ca.target_id = t.targetid
    JOIN targettobasetarget ttb ON t.targetid = ttb.target_id
    JOIN basetarget bt ON ttb.basetarget_id = bt.basetargetid
    LEFT JOIN actiontype at ON ca.actiontype_id = at.action
    WHERE c.name LIKE ? COLLATE NOCASE
    ORDER BY c.name, ca.primary_target DESC, ca.drug_target DESC
    """
    return run_query(sql, params=(name,))


def get_compound_sets(name):
    """
    Get all compound sets a compound belongs to, with set type labels.

    Path: compound → compoundtocompoundset → compoundset → compoundsettype

    Returns DataFrame with pdid, compound_name, set_name, set_type, version.
    """
    sql = """
    SELECT
        c.pdid,
        c.name AS compound_name,
        cs.name AS set_name,
        cst.label AS set_type,
        cs.version
    FROM compound c
    JOIN compoundtocompoundset ctcs ON c.compoundid = ctcs.compound_id
    JOIN compoundset cs ON ctcs.compoundset_id = cs.compoundsetid
    JOIN compoundsettype cst ON cs.compoundsettype_id = cst.compoundsettypeid
    WHERE c.name LIKE ? COLLATE NOCASE
    ORDER BY cst.label, cs.name
    """
    return run_query(sql, params=(name,))


def get_primary_target(name):
    """
    Get the primary target(s) for a compound (primary_target flag = 1).

    Returns DataFrame with pdid, compound_name, target_gene, action,
             primary_target, drug_target.
    """
    sql = """
    SELECT
        c.pdid,
        c.name AS compound_name,
        bt.gene_name AS target_gene,
        ca.actiontype_id AS action,
        ca.primary_target,
        ca.drug_target
    FROM compound c
    JOIN compoundaction ca ON c.compoundid = ca.compound_id
    JOIN target t ON ca.target_id = t.targetid
    JOIN targettobasetarget ttb ON t.targetid = ttb.target_id
    JOIN basetarget bt ON ttb.basetarget_id = bt.basetargetid
    WHERE c.name LIKE ? COLLATE NOCASE
      AND ca.primary_target = 1
    """
    return run_query(sql, params=(name,))


def get_selectivity(name):
    """
    Get pre-computed selectivity metrics from compoundbasetargetcriteria.

    Returns DataFrame with pdid, compound_name, target_gene, potency,
             selectivity, selectivity_score, family_selectivity, cell_potency.
    """
    sql = """
    SELECT
        c.pdid,
        c.name AS compound_name,
        bt.gene_name AS target_gene,
        cbtc.potency,
        cbtc.selectivity,
        cbtc.selectivity_score,
        cbtc.family_selectivity,
        cbtc.cell_potency,
        cbtc.potency_selectivity_synergy
    FROM compound c
    JOIN compoundbasetargetcriteria cbtc ON c.compoundid = cbtc.compound_id
    JOIN basetarget bt ON cbtc.basetarget_id = bt.basetargetid
    WHERE c.name LIKE ? COLLATE NOCASE
    ORDER BY cbtc.selectivity DESC
    """
    return run_query(sql, params=(name,))


def compute_s_score(name, threshold=6.0):
    """
    Compute a selectivity S-score from raw activity data.

    S-score = fraction of tested targets with best potency >= threshold.
    This is the Karaman et al. (2008) selectivity score: S(s) = n_hit / n_tested,
    where n_hit = targets with potency above the threshold and n_tested = total
    targets with activity data. Lower S-score = more selective compound.

    For example, S(1 µM) = S(threshold=6.0) counts the fraction of targets
    inhibited/bound at 1 µM or better. A perfectly selective compound that hits
    only its intended target scores 1/n_tested; a promiscuous compound approaches 1.0.

    Limitation: all tested targets are weighted equally in the denominator,
    regardless of how many targets were screened. A compound tested against only
    5 targets (all potent) gets S=1.0, while one tested against 300 targets with
    50 potent gets S=0.17 — the latter is actually more selective despite the
    higher score. Compare S-scores only between compounds with similar screening
    breadth, or use the Gini coefficient or partition index for cross-breadth
    comparisons.

    Parameters
    ----------
    name : str        — compound name (case-insensitive)
    threshold : float — potency cutoff (default 6.0 ≈ 1 µM)

    Returns dict with compound_name, n_targets, n_hit, s_score, top_target.
    """
    df = get_compound_targets(name, log_only=True)
    if df.empty:
        return {'compound_name': name, 'n_targets': 0, 'n_hit': 0, 's_score': 0, 'top_target': None}
    n_total = len(df)
    n_hit = int((df['best_potency'] >= threshold).sum())
    s_score = n_hit / n_total if n_total > 0 else 0
    top = df.iloc[0]
    return {
        'compound_name': name,
        'n_targets': n_total,
        'n_hit': n_hit,
        's_score': round(s_score, 3),
        'top_target': f"{top['target_gene']} ({top['best_activity_type']}={top['best_potency']})"
    }


# ════════════════════════════════════════════════════════════════
#  TARGET-CENTRIC QUERIES (US2)
# ════════════════════════════════════════════════════════════════

def get_basetarget(gene):
    """
    Look up basetarget(s) by gene name (case-insensitive).

    Returns DataFrame with basetargetid, gene_name, name, human_uniprot_id,
             target_family, targettype_name.
    """
    sql = """
    SELECT
        bt.basetargetid,
        bt.gene_name,
        bt.name,
        bt.human_uniprot_id,
        bt.target_family,
        tt.name AS targettype_name
    FROM basetarget bt
    LEFT JOIN targettype tt ON bt.targettype_id = tt.targettypeid
    WHERE bt.gene_name LIKE ? COLLATE NOCASE
    """
    return run_query(sql, params=(gene,))


def get_target_compounds(gene, inactive=False, log_only=True, min_confidence=1):
    """
    Get all compounds with activity against a target gene, with set info.

    Path: basetarget → targettobasetarget → target → activity → compound
          → compoundtocompoundset → compoundset → compoundsettype

    Parameters
    ----------
    gene : str            — target gene name
    inactive : bool       — if True, include inactive measurements
    log_only : bool       — if True (default), restrict to log-scale activity types
                            (pIC50, pKd, pKi, pEC50, pAC50, pPotency) to avoid mixing
                            with percentage-scale values (Inhibition 0-100%, Dmax)
    min_confidence : int  — minimum confidence level (1=direct, 2=derived). Default 1.

    Returns DataFrame with pdid, compound_name, n_sets, set_types.
    """
    inactive_filter = '' if inactive else 'AND a.inactive = 0'
    if log_only:
        quoted = ",".join(f'"{t}"' for t in LOG_ACTIVITY_TYPES)
        type_filter = f'AND a.activity_type IN ({quoted})'
    else:
        type_filter = ''
    conf_filter = f'AND a.confidence <= {int(min_confidence)}'
    sql = f"""
    SELECT
        c.pdid,
        c.name AS compound_name,
        COUNT(DISTINCT cs.compoundsetid) AS n_sets,
        GROUP_CONCAT(DISTINCT cst.label) AS set_types
    FROM basetarget bt
    JOIN targettobasetarget ttb ON bt.basetargetid = ttb.basetarget_id
    JOIN target t ON ttb.target_id = t.targetid
    JOIN activity a ON t.targetid = a.target_id
    JOIN compound c ON a.compound_id = c.compoundid
    LEFT JOIN compoundtocompoundset ctcs ON c.compoundid = ctcs.compound_id
    LEFT JOIN compoundset cs ON ctcs.compoundset_id = cs.compoundsetid
    LEFT JOIN compoundsettype cst ON cs.compoundsettype_id = cst.compoundsettypeid
    WHERE bt.gene_name LIKE ? COLLATE NOCASE
      {inactive_filter}
      {type_filter}
      {conf_filter}
    GROUP BY c.pdid, c.name
    ORDER BY n_sets DESC
    """
    return run_query(sql, params=(gene,))


def get_target_settype_breakdown(gene, log_only=True, min_confidence=1):
    """
    Get breakdown of compounds targeting a gene, grouped by set type.

    Parameters
    ----------
    gene : str            — target gene name
    log_only : bool       — if True (default), restrict to log-scale activity types
    min_confidence : int  — minimum confidence level (1=direct, 2=derived). Default 1.

    Returns DataFrame with set_type, n_compounds.
    """
    if log_only:
        quoted = ",".join(f'"{t}"' for t in LOG_ACTIVITY_TYPES)
        type_filter = f'AND a.activity_type IN ({quoted})'
    else:
        type_filter = ''
    conf_filter = f'AND a.confidence <= {int(min_confidence)}'
    sql = f"""
    SELECT
        cst.label AS set_type,
        COUNT(DISTINCT c.compoundid) AS n_compounds
    FROM basetarget bt
    JOIN targettobasetarget ttb ON bt.basetargetid = ttb.basetarget_id
    JOIN target t ON ttb.target_id = t.targetid
    JOIN activity a ON t.targetid = a.target_id
    JOIN compound c ON a.compound_id = c.compoundid
    JOIN compoundtocompoundset ctcs ON c.compoundid = ctcs.compound_id
    JOIN compoundset cs ON ctcs.compoundset_id = cs.compoundsetid
    JOIN compoundsettype cst ON cs.compoundsettype_id = cst.compoundsettypeid
    WHERE bt.gene_name LIKE ? COLLATE NOCASE
      AND a.inactive = 0
      {type_filter}
      {conf_filter}
    GROUP BY cst.label
    ORDER BY n_compounds DESC
    """
    return run_query(sql, params=(gene,))


def get_most_potent_compounds(gene, limit=10, log_only=True, exact_only=True,
                              min_confidence=1):
    """
    Rank compounds by best potency against a target gene.

    Parameters
    ----------
    gene : str            — target gene name
    limit : int           — number of top compounds to return
    log_only : bool       — if True (default), filter to log-scale activity types only
    exact_only : bool     — if True (default), exclude compounds whose ONLY measurements
                            use the '>' operator (screening negatives). These would
                            appear as "most potent" with misleadingly low values.
    min_confidence : int  — minimum confidence level (1=direct, 2=derived). Default 1.

    Returns DataFrame with pdid, compound_name, best_potency, activity_types, n_measurements.
    """
    if log_only:
        quoted_types = ",".join(f'"{t}"' for t in LOG_ACTIVITY_TYPES)
        type_filter = f"AND a.activity_type IN ({quoted_types})"
    else:
        type_filter = ''
    conf_filter = f'AND a.confidence <= {int(min_confidence)}'
    # When exact_only, restrict to '=' value_type (exclude '>' screening negatives)
    exact_filter = "AND a.value_type = '='" if exact_only else ''
    sql = f"""
    SELECT
        c.pdid,
        c.name AS compound_name,
        MAX(a.activity_value) AS best_potency,
        GROUP_CONCAT(DISTINCT a.activity_type) AS activity_types,
        COUNT(*) AS n_measurements
    FROM basetarget bt
    JOIN targettobasetarget ttb ON bt.basetargetid = ttb.basetarget_id
    JOIN target t ON ttb.target_id = t.targetid
    JOIN activity a ON t.targetid = a.target_id
    JOIN compound c ON a.compound_id = c.compoundid
    WHERE bt.gene_name LIKE ? COLLATE NOCASE
      AND a.inactive = 0
      {conf_filter}
      AND a.activity_value IS NOT NULL
      {type_filter}
      {exact_filter}
    GROUP BY c.pdid, c.name
    ORDER BY best_potency DESC
    LIMIT ?
    """
    return run_query(sql, params=(gene, limit))


def get_most_selective_compounds(gene, limit=10):
    """
    Rank compounds by pre-computed selectivity against a target gene.

    Returns DataFrame with pdid, compound_name, potency, selectivity,
             selectivity_score, family_selectivity.
    """
    sql = """
    SELECT
        c.pdid,
        c.name AS compound_name,
        cbtc.potency,
        cbtc.selectivity,
        cbtc.selectivity_score,
        cbtc.family_selectivity,
        cbtc.potency_selectivity_synergy
    FROM basetarget bt
    JOIN compoundbasetargetcriteria cbtc ON bt.basetargetid = cbtc.basetarget_id
    JOIN compound c ON cbtc.compound_id = c.compoundid
    WHERE bt.gene_name LIKE ? COLLATE NOCASE
      AND cbtc.selectivity IS NOT NULL
    ORDER BY cbtc.selectivity DESC
    LIMIT ?
    """
    return run_query(sql, params=(gene, limit))


def get_potency_selectivity(gene):
    """
    Get potency and selectivity pairs for all compounds targeting a gene.

    Returns DataFrame with pdid, compound_name, potency, selectivity,
             potency_selectivity_synergy.
    """
    sql = """
    SELECT
        c.pdid,
        c.name AS compound_name,
        cbtc.potency,
        cbtc.selectivity,
        cbtc.potency_selectivity_synergy
    FROM basetarget bt
    JOIN compoundbasetargetcriteria cbtc ON bt.basetargetid = cbtc.basetarget_id
    JOIN compound c ON cbtc.compound_id = c.compoundid
    WHERE bt.gene_name LIKE ? COLLATE NOCASE
      AND cbtc.potency IS NOT NULL
      AND cbtc.selectivity IS NOT NULL
    """
    return run_query(sql, params=(gene,), show_raw=False)


def get_target_probes(gene):
    """
    Get probes targeting a gene, with their control compounds.

    Uses the proper probetobasetarget + probecontrol tables.

    Path: basetarget → probetobasetarget → probe → compound
          probe → probecontrol → compound (control)

    Returns DataFrame with target_gene, probe_pdid, probe_name, probe_origin,
             control_name, control_compound_id, obsolete_flag.
    """
    sql = """
    SELECT DISTINCT
        bt.gene_name AS target_gene,
        cp.pdid AS probe_pdid,
        cp.name AS probe_name,
        p.origin_id AS probe_origin,
        p.obsolete_flag,
        p.control_flag,
        pc.name AS control_name,
        pc.compound_id AS control_compound_id,
        pc.smiles AS control_smiles
    FROM basetarget bt
    JOIN probetobasetarget pbt ON bt.basetargetid = pbt.basetarget_id
    JOIN probe p ON pbt.probe_id = p.probeid
    JOIN compound cp ON p.compound_id = cp.compoundid
    LEFT JOIN probecontrol pc ON pc.probe_id = p.probeid
    WHERE bt.gene_name LIKE ? COLLATE NOCASE
    ORDER BY p.obsolete_flag, cp.name, pc.name
    """
    return run_query(sql, params=(gene,))


def get_murcko_scaffolds(gene=None, compound_ids=None, limit=100, min_confidence=1):
    """
    Get pre-computed Murcko scaffolds for compounds targeting a gene.

    Uses the compoundtoscaffold + scaffold + scaffoldtype tables
    (pre-computed, no RDKit needed).

    Parameters
    ----------
    gene : str            — target gene name (alternative to compound_ids)
    compound_ids : list   — list of compound IDs (alternative to gene)
    limit : int           — max compounds to analyze
    min_confidence : int  — minimum confidence level (1=direct, 2=derived). Default 1.
    """
    if gene is not None:
        sql = f"""
        SELECT
            c.pdid,
            c.name AS compound_name,
            MAX(a.activity_value) AS best_potency,
            s.smiles AS scaffold_smiles
        FROM basetarget bt
        JOIN targettobasetarget ttb ON bt.basetargetid = ttb.basetarget_id
        JOIN target t ON ttb.target_id = t.targetid
        JOIN activity a ON t.targetid = a.target_id
        JOIN compound c ON a.compound_id = c.compoundid
        JOIN compoundtoscaffold cts ON c.compoundid = cts.compound_id
        JOIN scaffold s ON cts.scaffold_id = s.scaffoldid
        JOIN scaffoldtype st ON s.scaffoldtype_id = st.scaffoldtypeid
        WHERE bt.gene_name LIKE ? COLLATE NOCASE
          AND a.inactive = 0
          AND a.confidence <= {int(min_confidence)}
          AND a.activity_value IS NOT NULL
          AND a.activity_type IN ({','.join(f'"{t}"' for t in LOG_ACTIVITY_TYPES)})
          AND st.name = 'Murcko'
          AND s.smiles != ''
        GROUP BY c.pdid, c.name, s.smiles
        ORDER BY best_potency DESC
        LIMIT ?
        """
        return run_query(sql, params=(gene, limit), show_raw=False)
    elif compound_ids is not None:
        placeholders = ','.join('?' * len(compound_ids))
        sql = f"""
        SELECT
            c.pdid,
            c.name AS compound_name,
            s.smiles AS scaffold_smiles
        FROM compound c
        JOIN compoundtoscaffold cts ON c.compoundid = cts.compound_id
        JOIN scaffold s ON cts.scaffold_id = s.scaffoldid
        JOIN scaffoldtype st ON s.scaffoldtype_id = st.scaffoldtypeid
        WHERE c.compoundid IN ({placeholders})
          AND st.name = 'Murcko'
          AND s.smiles != ''
        """
        return run_query(sql, params=compound_ids, show_raw=False)
    else:
        raise ValueError("Provide either gene or compound_ids")


def get_scaffold_frequency(scaffold_df):
    """
    Count scaffold frequency from a scaffold DataFrame.

    Returns DataFrame sorted by frequency: scaffold_smiles, n_compounds.
    """
    counts = Counter(scaffold_df['scaffold_smiles'].dropna())
    return pd.DataFrame(
        counts.most_common(),
        columns=['scaffold_smiles', 'n_compounds']
    )


# ════════════════════════════════════════════════════════════════
#  PROTEIN FAMILY QUERIES (US3)
# ════════════════════════════════════════════════════════════════

def get_family_members(gene_names):
    """
    Find all basetargets whose gene_name includes any of the given genes.

    Captures both single-gene targets (e.g. 'KRAS') and composite targets
    (e.g. 'KRAS,SOS1', 'HRAS,KRAS,NRAS').

    Parameters
    ----------
    gene_names : list — e.g. ['KRAS', 'NRAS', 'HRAS']

    Returns DataFrame with basetarget_id, gene_name, basetarget_name,
             uniprot_id, target_family, targettype_name, target_type.
    """
    # Build WHERE clause: exact match OR comma-delimited match
    conditions = []
    params = []
    for g in gene_names:
        conditions.append(f"bt.gene_name = ?")
        conditions.append(f"bt.gene_name LIKE ?")
        conditions.append(f"bt.gene_name LIKE ?")
        conditions.append(f"bt.gene_name LIKE ?")
        params.extend([g, f'{g},%', f'%,{g}', f'%,{g},%'])

    where = ' OR '.join(conditions)
    sql = f"""
    SELECT DISTINCT
        bt.basetargetid AS basetarget_id,
        bt.gene_name,
        bt.name AS basetarget_name,
        bt.human_uniprot_id AS uniprot_id,
        bt.target_family,
        tt.name AS targettype_name
    FROM basetarget bt
    LEFT JOIN targettype tt ON bt.targettype_id = tt.targettypeid
    WHERE {where}
    ORDER BY bt.gene_name
    """
    df = run_query(sql, params=params, show_raw=False)
    df['target_type'] = df['gene_name'].apply(
        lambda g: 'single-gene' if g in gene_names else 'composite'
    )
    return df


def get_family_compounds(gene_names, inactive=False, log_only=True, min_confidence=1):
    """
    Get all compounds with activity against any member of a protein family.

    Parameters
    ----------
    gene_names : list     — e.g. ['KRAS', 'NRAS', 'HRAS']
    inactive : bool       — if True, include inactive measurements
    log_only : bool       — if True (default), restrict to log-scale activity types
                            (pIC50, pKd, pKi, pEC50, pAC50, pPotency) to avoid mixing
                            with percentage-scale values (Inhibition 0-100%, Dmax)
    min_confidence : int  — minimum confidence level (1=direct, 2=derived). Default 1.

    Returns DataFrame with target_gene, pdid, compound_name, activity_type,
             activity_value, set_name, set_type_label.
    """
    members = get_family_members(gene_names)
    if members.empty:
        return pd.DataFrame()
    bt_ids = members['basetarget_id'].tolist()
    placeholders = ','.join('?' * len(bt_ids))
    inactive_filter = '' if inactive else 'AND a.inactive = 0'
    if log_only:
        quoted = ",".join(f'"{t}"' for t in LOG_ACTIVITY_TYPES)
        type_filter = f'AND a.activity_type IN ({quoted})'
    else:
        type_filter = ''
    conf_filter = f'AND a.confidence <= {int(min_confidence)}'

    sql = f"""
    SELECT DISTINCT
        bt.gene_name AS target_gene,
        c.pdid,
        c.name AS compound_name,
        a.activity_type,
        a.activity_value,
        cs.name AS set_name,
        cst.label AS set_type_label
    FROM basetarget bt
    JOIN targettobasetarget ttbt ON ttbt.basetarget_id = bt.basetargetid
    JOIN target t ON t.targetid = ttbt.target_id
    JOIN activity a ON a.target_id = t.targetid
    JOIN compound c ON c.compoundid = a.compound_id
    LEFT JOIN compoundtocompoundset ctcs ON ctcs.compound_id = c.compoundid
    LEFT JOIN compoundset cs ON cs.compoundsetid = ctcs.compoundset_id
    LEFT JOIN compoundsettype cst ON cst.compoundsettypeid = cs.compoundsettype_id
    WHERE bt.basetargetid IN ({placeholders})
      {inactive_filter}
      {type_filter}
      {conf_filter}
    ORDER BY bt.gene_name, c.name
    """
    return run_query(sql, params=bt_ids, show_raw=False)


def classify_compound_category(set_type_label, set_name):
    """
    Classify a compound-set membership into Probe, Drug, Chemogenomic, or Other.

    Parameters
    ----------
    set_type_label : str — from compoundsettype.label
    set_name : str       — from compoundset.name

    Returns str: 'Probe', 'Drug', 'Chemogenomic', or 'Other'
    """
    label = str(set_type_label or '')
    name = str(set_name or '')
    if any(kw.lower() in name.lower() for kw in CHEMOGEN_KEYWORDS):
        return 'Chemogenomic'
    if label == 'Probe compound sets':
        return 'Probe'
    if label == 'Drug compound sets':
        return 'Drug'
    return 'Other'


def get_family_compound_categories(gene_names):
    """
    Get per-compound category classification for a protein family.

    Returns DataFrame with pdid, compound_name, Probe, Drug, Chemogenomic, Other
    (binary flags — a compound can be in multiple categories).
    """
    df = get_family_compounds(gene_names)
    if df.empty:
        return pd.DataFrame()
    df_sets = df.drop_duplicates(subset=['pdid', 'compound_name', 'set_name']).copy()
    df_sets['category'] = df_sets.apply(
        lambda r: classify_compound_category(r['set_type_label'], r['set_name']), axis=1
    )
    cats = df_sets.groupby(['pdid', 'compound_name', 'category']).size().reset_index(name='n')
    pivot = cats.pivot_table(index=['pdid', 'compound_name'], columns='category', values='n', fill_value=0)
    for cat in ['Probe', 'Drug', 'Chemogenomic', 'Other']:
        if cat not in pivot.columns:
            pivot[cat] = 0
    return (pivot > 0).astype(int).reset_index()


def get_family_probes(gene_names):
    """
    Get probes targeting any member of a protein family, with controls.

    Uses probetobasetarget + probecontrol tables.

    Returns DataFrame with target_gene, probe_pdid, probe_name, probe_origin,
             control_name, obsolete_flag.
    """
    members = get_family_members(gene_names)
    if members.empty:
        return pd.DataFrame()
    bt_ids = members['basetarget_id'].tolist()
    placeholders = ','.join('?' * len(bt_ids))

    sql = f"""
    SELECT DISTINCT
        bt.gene_name AS target_gene,
        cp.pdid AS probe_pdid,
        cp.name AS probe_name,
        p.origin_id AS probe_origin,
        p.obsolete_flag,
        pc.name AS control_name,
        pc.compound_id AS control_compound_id
    FROM basetarget bt
    JOIN probetobasetarget pbt ON bt.basetargetid = pbt.basetarget_id
    JOIN probe p ON pbt.probe_id = p.probeid
    JOIN compound cp ON p.compound_id = cp.compoundid
    LEFT JOIN probecontrol pc ON pc.probe_id = p.probeid
    WHERE bt.basetargetid IN ({placeholders})
    ORDER BY bt.gene_name, p.obsolete_flag, cp.name
    """
    return run_query(sql, params=bt_ids)


# ════════════════════════════════════════════════════════════════
#  TARGET COVERAGE QUERIES (US4)
# ════════════════════════════════════════════════════════════════

def get_target_coverage_top(limit=20, inactive=False, log_only=True, min_confidence=1):
    """
    Get the most promiscuous compounds (most distinct basetargets).

    Parameters
    ----------
    limit : int           — number of top compounds to return
    inactive : bool       — if True, include inactive measurements
    log_only : bool       — if True (default), restrict to log-scale activity types
                            to avoid counting percentage-scale measurements
    min_confidence : int  — minimum confidence level (1=direct, 2=derived). Default 1.

    Returns DataFrame with pdid, compound_name, n_targets.
    """
    inactive_filter = '' if inactive else 'AND a.inactive = 0'
    if log_only:
        quoted = ",".join(f'"{t}"' for t in LOG_ACTIVITY_TYPES)
        type_filter = f'AND a.activity_type IN ({quoted})'
    else:
        type_filter = ''
    conf_filter = f'AND a.confidence <= {int(min_confidence)}'
    sql = f"""
    SELECT
        c.pdid,
        c.name AS compound_name,
        COUNT(DISTINCT bt.basetargetid) AS n_targets
    FROM compound c
    JOIN activity a ON a.compound_id = c.compoundid
    JOIN target t ON t.targetid = a.target_id
    JOIN targettobasetarget ttbt ON ttbt.target_id = t.targetid
    JOIN basetarget bt ON bt.basetargetid = ttbt.basetarget_id
    WHERE 1=1
      {inactive_filter}
      {type_filter}
      {conf_filter}
    GROUP BY c.pdid, c.name
    ORDER BY n_targets DESC
    LIMIT ?
    """
    return run_query(sql, params=(limit,))


def get_set_coverage(set_name, inactive=False, log_only=True, min_confidence=1):
    """
    Get target coverage for a single compound set.

    Parameters
    ----------
    set_name : str        — compound set name
    inactive : bool       — if True, include inactive measurements
    log_only : bool       — if True (default), restrict to log-scale activity types
    min_confidence : int  — minimum confidence level (1=direct, 2=derived). Default 1.

    Returns DataFrame with gene_name, target_name, target_family, n_compounds.
    """
    inactive_filter = '' if inactive else 'AND a.inactive = 0'
    if log_only:
        quoted = ",".join(f'"{t}"' for t in LOG_ACTIVITY_TYPES)
        type_filter = f'AND a.activity_type IN ({quoted})'
    else:
        type_filter = ''
    conf_filter = f'AND a.confidence <= {int(min_confidence)}'
    sql = f"""
    SELECT
        bt.gene_name,
        bt.name AS target_name,
        bt.target_family,
        COUNT(DISTINCT c.pdid) AS n_compounds
    FROM compoundset cs
    JOIN compoundtocompoundset ctcs ON ctcs.compoundset_id = cs.compoundsetid
    JOIN compound c ON c.compoundid = ctcs.compound_id
    JOIN activity a ON a.compound_id = c.compoundid
    JOIN target t ON t.targetid = a.target_id
    JOIN targettobasetarget ttbt ON ttbt.target_id = t.targetid
    JOIN basetarget bt ON bt.basetargetid = ttbt.basetarget_id
    WHERE cs.name = ?
      {inactive_filter}
      {type_filter}
      {conf_filter}
    GROUP BY bt.gene_name, bt.name, bt.target_family
    ORDER BY n_compounds DESC
    """
    return run_query(sql, params=(set_name,))


def get_set_compounds(set_name):
    """
    Get all compounds in a named compound set.

    Returns DataFrame with pdid, compound_name, set_type.
    """
    sql = """
    SELECT
        c.pdid,
        c.name AS compound_name,
        cst.label AS set_type
    FROM compoundset cs
    JOIN compoundtocompoundset ctcs ON ctcs.compoundset_id = cs.compoundsetid
    JOIN compound c ON c.compoundid = ctcs.compound_id
    JOIN compoundsettype cst ON cs.compoundsettype_id = cst.compoundsettypeid
    WHERE cs.name = ?
    ORDER BY c.name
    """
    return run_query(sql, params=(set_name,), show_raw=False)


def get_set_intersection(set_names):
    """
    Get compounds that appear in ALL of the specified sets.

    Parameters
    ----------
    set_names : list — e.g. ['High-quality chemical probes', 'SGC Probes']

    Returns DataFrame with pdid, compound_name, set_names, n_sets.
    """
    placeholders = ','.join('?' * len(set_names))
    sql = f"""
    SELECT
        c.pdid,
        c.name AS compound_name,
        GROUP_CONCAT(DISTINCT cs.name) AS set_names,
        COUNT(DISTINCT cs.compoundsetid) AS n_sets
    FROM compound c
    JOIN compoundtocompoundset ctcs ON ctcs.compound_id = c.compoundid
    JOIN compoundset cs ON cs.compoundsetid = ctcs.compoundset_id
    WHERE cs.name IN ({placeholders})
    GROUP BY c.pdid, c.name
    HAVING COUNT(DISTINCT cs.compoundsetid) >= ?
    ORDER BY c.name
    """
    return run_query(sql, params=list(set_names) + [len(set_names)])


def get_set_sizes(set_names):
    """
    Get compound counts for each set individually.

    Returns DataFrame with name, n_compounds.
    """
    placeholders = ','.join('?' * len(set_names))
    sql = f"""
    SELECT cs.name, COUNT(DISTINCT ctcs.compound_id) AS n_compounds
    FROM compoundset cs
    JOIN compoundtocompoundset ctcs ON ctcs.compoundset_id = cs.compoundsetid
    WHERE cs.name IN ({placeholders})
    GROUP BY cs.name
    """
    return run_query(sql, params=set_names)


def get_compound_scores(gene, score_name=None):
    """
    Get pre-computed compound-target scores from the score table.

    Parameters
    ----------
    gene : str            — target gene name
    score_name : str|None — filter to a specific score (e.g. 'Probe Miner Score').
                            If None, returns all scores.

    Returns DataFrame with pdid, compound_name, score_name, value, percentage.
    """
    score_filter = "AND s.name = ?" if score_name else ""
    params = [gene] + ([score_name] if score_name else [])
    sql = f"""
    SELECT
        c.pdid,
        c.name AS compound_name,
        s.name AS score_name,
        s.acronym AS score_acronym,
        cts.value,
        cts.percentage
    FROM basetarget bt
    JOIN compoundtargetscore cts ON bt.basetargetid = cts.basetarget_id
    JOIN score s ON cts.score_id = s.scoreid
    JOIN compound c ON cts.compound_id = c.compoundid
    WHERE bt.gene_name LIKE ? COLLATE NOCASE
      {score_filter}
    ORDER BY cts.value DESC
    """
    return run_query(sql, params=params)


def list_score_types():
    """List all available score types in the database."""
    sql = "SELECT scoreid, name, acronym, description FROM score ORDER BY scoreid"
    return run_query(sql, show_raw=False)


def list_compoundset_types():
    """List all compound set types with compound counts."""
    sql = """
    SELECT cst.label AS set_type, COUNT(DISTINCT ctcs.compound_id) AS n_compounds
    FROM compoundsettype cst
    JOIN compoundset cs ON cs.compoundsettype_id = cst.compoundsettypeid
    JOIN compoundtocompoundset ctcs ON ctcs.compoundset_id = cs.compoundsetid
    GROUP BY cst.label
    ORDER BY n_compounds DESC
    """
    return run_query(sql, show_raw=False)


# ════════════════════════════════════════════════════════════════
#  PLOTTING HELPERS
# ════════════════════════════════════════════════════════════════

def plot_potency_profile(df, compound_name, primary_genes=None, threshold=6.0,
                         save_path=None, hue_by_type=False, detailed_df=None):
    """
    Horizontal bar chart of compound potency across all targets.

    Parameters
    ----------
    df : DataFrame       — from get_compound_targets() (must have target_gene,
                           best_potency, best_activity_type)
    compound_name : str  — for title
    primary_genes : set  — genes to highlight with a marker
    threshold : float    — vertical line at this potency value
    save_path : str|None — if given, save figure to this path
    hue_by_type : bool   — if True, color bars by activity_type instead of
                           primary/off-target. Requires detailed_df.
    detailed_df : DataFrame — from get_compound_targets_detailed(), used when
                              hue_by_type=True to draw grouped bars per target.
    """
    primary_genes = primary_genes or set()

    if hue_by_type and detailed_df is not None:
        _plot_potency_profile_grouped(detailed_df, compound_name, threshold,
                                      save_path)
        return

    df = df.sort_values('best_potency', ascending=True).copy()

    # Color by primary vs off-target
    colors = [PALETTE[0] if g in primary_genes else PALETTE[5]
              for g in df['target_gene']]

    fig, ax = plt.subplots(figsize=(11, max(4, len(df) * 0.4)))
    bars = ax.barh(df['target_gene'], df['best_potency'], color=colors,
                   edgecolor='white', height=0.7)
    ax.axvline(x=threshold, color=PALETTE[1], linestyle='--', linewidth=1.2)

    # Apply hatching to contradictory bars and annotate
    has_flag_col = 'contradiction_flag' in df.columns
    for bar, (_, row) in zip(bars, df.iterrows()):
        atype = row.get('best_activity_type', '?')
        val = row['best_potency']
        label = f'{atype}={val:.1f}'
        if has_flag_col and row.get('contradiction_flag', False):
            bar.set_hatch('///')
            bar.set_edgecolor(PALETTE[1])
            label += ' (!)'
        ax.text(val + 0.08, bar.get_y() + bar.get_height() / 2,
                label, va='center', fontsize=8, color='#333333')

    ax.set_xlabel('Best Potency (-log$_{10}$ M)')
    ax.set_ylabel('Target')
    ax.set_title(f'{compound_name}: Potency Profile Across Targets')

    from matplotlib.patches import Patch
    legend = [Patch(facecolor=PALETTE[0], label='Primary targets'),
              Patch(facecolor=PALETTE[5], label='Off-targets'),
              plt.Line2D([0], [0], color=PALETTE[1], linestyle='--',
                         label=f'p = {threshold} ({10**-threshold * 1e6:.0f} µM)')]
    if has_flag_col and df['contradiction_flag'].any():
        legend.append(Patch(facecolor=PALETTE[5], hatch='///',
                            edgecolor=PALETTE[1],
                            label='Contradictory data (!)'))
    ax.legend(handles=legend, loc='lower right', fontsize=9)
    plt.tight_layout()
    if save_path:
        plt.savefig(save_path, dpi=150, bbox_inches='tight')
    plt.show()


def _plot_potency_profile_grouped(detailed_df, compound_name, threshold=6.0,
                                  save_path=None):
    """Internal: grouped horizontal bar chart with one bar per activity type per target."""
    # Pivot: rows=target_gene, columns=activity_type, values=best_potency
    pivot = detailed_df.pivot_table(index='target_gene', columns='activity_type',
                                    values='best_potency', aggfunc='max')
    # Sort targets by their overall max potency (ascending so strongest is at top)
    max_per_target = pivot.max(axis=1).sort_values(ascending=True)
    pivot = pivot.loc[max_per_target.index]

    # Assign colors per activity type
    activity_types = list(pivot.columns)
    type_colors = {at: PALETTE[i % len(PALETTE)] for i, at in enumerate(activity_types)}

    n_targets = len(pivot)
    n_types = len(activity_types)
    bar_height = 0.8 / n_types

    fig, ax = plt.subplots(figsize=(11, max(5, n_targets * 0.5)))

    for i, atype in enumerate(activity_types):
        offsets = np.arange(n_targets) + i * bar_height - 0.4 + bar_height / 2
        vals = pivot[atype].values
        # NaN → 0 for plotting, but we'll mask them
        masked_vals = np.nan_to_num(vals, nan=0)
        bars = ax.barh(offsets, masked_vals, height=bar_height,
                       color=type_colors[atype], edgecolor='white', linewidth=0.4,
                       label=atype)
        # Annotate non-zero bars with value
        for bar, val in zip(bars, vals):
            if not np.isnan(val) and val > 0:
                ax.text(val + 0.08, bar.get_y() + bar.get_height() / 2,
                        f'{val:.1f}', va='center', fontsize=7, color='#555555')

    ax.set_yticks(np.arange(n_targets))
    ax.set_yticklabels(pivot.index)
    ax.axvline(x=threshold, color=PALETTE[1], linestyle='--', linewidth=1.2)
    ax.set_xlabel('Best Potency (-log$_{10}$ M)')
    ax.set_ylabel('Target')
    ax.set_title(f'{compound_name}: Potency by Measurement Type')

    from matplotlib.patches import Patch
    type_handles = [Patch(facecolor=type_colors[at], label=at)
                    for at in activity_types]
    type_handles.append(
        plt.Line2D([0], [0], color=PALETTE[1], linestyle='--',
                   label=f'p = {threshold} ({10**-threshold * 1e6:.0f} µM)'))
    ax.legend(handles=type_handles, loc='lower right', fontsize=8,
              title='Activity type', title_fontsize=9)
    plt.tight_layout()
    if save_path:
        plt.savefig(save_path, dpi=150, bbox_inches='tight')
    plt.show()


def plot_potency_vs_selectivity(df, title='Potency vs Selectivity', save_path=None):
    """
    Scatter plot of potency vs selectivity for compounds targeting a gene.

    Parameters
    ----------
    df : DataFrame — from get_potency_selectivity()
    title : str
    save_path : str|None
    """
    fig, ax = plt.subplots(figsize=(9, 6))
    colors = [PALETTE[2] if s == 1 else PALETTE[5] for s in df['potency_selectivity_synergy']]
    ax.scatter(df['potency'], df['selectivity'], c=colors, alpha=0.7, s=50,
               edgecolors='white', linewidth=0.5)

    top = df.nlargest(5, 'selectivity')
    for _, row in top.iterrows():
        ax.annotate(str(row['compound_name'])[:15], (row['potency'], row['selectivity']),
                    fontsize=7, ha='left', va='bottom', xytext=(3, 3), textcoords='offset points')

    ax.set_xlabel('Curated Potency (-log$_{10}$ M)')
    ax.set_ylabel('Selectivity Window (log units)')
    ax.set_title(title)

    from matplotlib.patches import Patch
    ax.legend(handles=[
        Patch(facecolor=PALETTE[2], label='Potency-selectivity synergy'),
        Patch(facecolor=PALETTE[5], label='No synergy'),
    ], loc='upper left', fontsize=9)

    plt.tight_layout()
    if save_path:
        plt.savefig(save_path, dpi=150, bbox_inches='tight')
    plt.show()


def plot_scaffold_frequency(scaffold_counts_df, title='Top Murcko Scaffolds',
                            top_n=10, save_path=None):
    """
    Horizontal bar chart of scaffold frequency.

    Parameters
    ----------
    scaffold_counts_df : DataFrame — from get_scaffold_frequency()
    title : str
    top_n : int
    save_path : str|None
    """
    top = scaffold_counts_df.head(top_n)
    labels = [s[:40] + '...' if len(s) > 40 else s for s in top['scaffold_smiles']]

    fig, ax = plt.subplots(figsize=(10, 5))
    ax.barh(range(len(labels)), top['n_compounds'], color=PALETTE[0],
            edgecolor='white', height=0.7)
    ax.set_yticks(range(len(labels)))
    ax.set_yticklabels(labels, fontsize=8)
    ax.invert_yaxis()
    ax.set_xlabel('Number of Compounds')
    ax.set_title(title)
    for i, count in enumerate(top['n_compounds']):
        ax.text(count + 0.3, i, str(count), va='center', fontsize=9)
    plt.tight_layout()
    if save_path:
        plt.savefig(save_path, dpi=150, bbox_inches='tight')
    plt.show()


def plot_stacked_categories(df, title='Compounds by Target and Category',
                            save_path=None):
    """
    Stacked horizontal bar chart of compound categories per target gene.

    Parameters
    ----------
    df : DataFrame — from get_family_compound_categories() joined with target_gene
    title : str
    save_path : str|None
    """
    categories = ['Probe', 'Drug', 'Chemogenomic', 'Other']
    colors = [PALETTE[0], PALETTE[1], PALETTE[2], PALETTE[5]]

    # Sum per target_gene
    if 'target_gene' in df.columns:
        summary = df.groupby('target_gene')[categories].sum()
    else:
        summary = df[categories].to_frame().T

    summary = summary.sort_values('Other', ascending=False)

    fig, ax = plt.subplots(figsize=(10, max(4, len(summary) * 0.5)))
    bottom = np.zeros(len(summary))
    for cat, color in zip(categories, colors):
        vals = summary[cat].values
        ax.barh(summary.index, vals, left=bottom, label=cat, color=color,
                edgecolor='white', linewidth=0.5)
        bottom += vals

    ax.set_xlabel('Number of Distinct Compounds')
    ax.set_title(title)
    ax.legend(loc='lower right', frameon=True)
    ax.invert_yaxis()
    plt.tight_layout()
    if save_path:
        plt.savefig(save_path, dpi=150, bbox_inches='tight')
    plt.show()


def plot_set_intersection(set_sizes_df, intersection_df, title='Probe Set Sizes and Intersection',
                          save_path=None):
    """
    Bar chart of set sizes, intersection, and union.

    Parameters
    ----------
    set_sizes_df : DataFrame      — from get_set_sizes()
    intersection_df : DataFrame   — from get_set_intersection()
    title : str
    save_path : str|None
    """
    fig, ax = plt.subplots(figsize=(8, 5))
    set_names = set_sizes_df['name'].tolist()
    set_vals = set_sizes_df['n_compounds'].tolist()
    intersection_n = len(intersection_df)
    union_n = sum(set_vals) - intersection_n

    categories = [s.replace(' ', '\n', 1) for s in set_names] + ['Intersection', 'Union']
    values = set_vals + [intersection_n, union_n]
    colors = [PALETTE[0], PALETTE[1], PALETTE[2], PALETTE[5]][:len(values)]

    bars = ax.bar(categories, values, color=colors, edgecolor='white', linewidth=0.5)
    ax.set_ylabel('Number of Compounds')
    ax.set_title(title)
    for bar, val in zip(bars, values):
        ax.text(bar.get_x() + bar.get_width() / 2, bar.get_height() + 10,
                str(val), ha='center', va='bottom', fontsize=11, fontweight='bold')
    plt.tight_layout()
    if save_path:
        plt.savefig(save_path, dpi=150, bbox_inches='tight')
    plt.show()


# ── Convenience: import all public names ─────────────────────────
__all__ = [
    # Core
    'run_query', 'get_connection', 'list_tables',
    # Compound-centric
    'get_compound', 'get_compound_targets', 'get_compound_targets_detailed',
    'get_compound_actions',
    'get_compound_sets', 'get_primary_target', 'get_selectivity', 'compute_s_score',
    # Target-centric
    'get_basetarget', 'get_target_compounds', 'get_target_settype_breakdown',
    'get_most_potent_compounds', 'get_most_selective_compounds',
    'get_potency_selectivity', 'get_target_probes',
    'get_murcko_scaffolds', 'get_scaffold_frequency',
    # Family-centric
    'get_family_members', 'get_family_compounds',
    'classify_compound_category', 'get_family_compound_categories',
    'get_family_probes',
    # Coverage
    'get_target_coverage_top', 'get_set_coverage', 'get_set_compounds',
    'get_set_intersection', 'get_set_sizes',
    'get_compound_scores', 'list_score_types', 'list_compoundset_types',
    # Plotting
    'plot_potency_profile', 'plot_potency_vs_selectivity',
    'plot_scaffold_frequency', 'plot_stacked_categories', 'plot_set_intersection',
    # Constants
    'PALETTE', 'LOG_ACTIVITY_TYPES', 'CHEMOGEN_KEYWORDS', 'DB_PATH',
]
