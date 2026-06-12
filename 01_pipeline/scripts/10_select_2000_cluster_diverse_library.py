from __future__ import annotations

from pathlib import Path
import argparse
import re
from collections import defaultdict
import numpy as np
import pandas as pd

BASE = Path.cwd()
SEQ = "utr5_sequence_tss_corrected"
DEFAULT_INPUTS = [
    BASE / "04_te_labeling/tables/tss_corrected_5utr_with_seq_clusters_and_heavy_scores.csv",
    BASE / "04_te_labeling/tables/tss_corrected_5utr_with_seq_clusters.csv",
    BASE / "04_te_labeling/tables/tss_corrected_5utr_multiomics_labels.csv",
    BASE / "04_te_labeling/tables/tss_corrected_5utr_robust_public_te_labels.csv",
]
OUT_CSV = BASE / "07_library_design/tables/selected_2000_50_100bp_cluster_diverse_evidence_balanced_library.csv"
OUT_FASTA = BASE / "07_library_design/fasta/selected_2000_50_100bp_cluster_diverse_evidence_balanced_library.fasta"
OUT_QC = BASE / "07_library_design/qc/selected_2000_50_100bp_cluster_diverse_evidence_balanced_summary.txt"
OUT_DIVERSITY = BASE / "06_modeling/tables/final_library_gene_cluster_diversity_summary.txt"
OUT_UAUG_SUMMARY_TXT = BASE / "07_library_design/qc/uaug_source_by_group_summary.txt"
OUT_UAUG_SUMMARY_CSV = BASE / "07_library_design/tables/uaug_source_by_group_summary.csv"
OUT_UAUG_POSITIVE_CSV = BASE / "07_library_design/tables/uaug_positive_final_library_rows.csv"
OUT_UAUG0_DRY_SUMMARY = BASE / "07_library_design/tables/uaug0_hard_filter_dry_run_summary.csv"
OUT_UAUG0_SHORTFALL = BASE / "07_library_design/tables/uaug0_hard_filter_quota_shortfall.csv"
OUT_UAUG0_REPLACEMENTS = BASE / "07_library_design/tables/uaug0_replacement_candidates.csv"
OUT_UAUG0_VALIDATION = BASE / "07_library_design/qc/uaug0_production_validation_report.txt"
OUT_REFILL_AUDIT = BASE / "07_library_design/tables/evidence_refill_audit.csv"
OUT_SELECTION_QC = BASE / "07_library_design/tables/v1.4_selection_policy_qc.csv"
for p in [OUT_CSV.parent, OUT_FASTA.parent, OUT_QC.parent, OUT_DIVERSITY.parent, OUT_UAUG_SUMMARY_TXT.parent, OUT_UAUG_SUMMARY_CSV.parent]:
    p.mkdir(parents=True, exist_ok=True)


def clean_seq(x):
    if pd.isna(x): return ""
    return re.sub(r"[^ACGTN]", "", str(x).upper().replace("U", "T"))


def gc_content(seq):
    seq = clean_seq(seq)
    return (seq.count("G") + seq.count("C")) / len(seq) if seq else np.nan


def choose_input(path_arg=None):
    if path_arg:
        p = Path(path_arg)
        if p.exists(): return p
        raise SystemExit(f"Input not found: {p}")
    for p in DEFAULT_INPUTS:
        if p.exists(): return p
    raise SystemExit("No label table found")


def forbidden_sites(seq):
    sites = {
        "BsaI_GGTCTC": "GGTCTC", "BsaI_GAGACC": "GAGACC",
        "BsmBI_CGTCTC": "CGTCTC", "BsmBI_GAGACG": "GAGACG",
        "EcoRI_GAATTC": "GAATTC", "XhoI_CTCGAG": "CTCGAG",
        "NheI_GCTAGC": "GCTAGC", "AgeI_ACCGGT": "ACCGGT", "NotI_GCGGCCGC": "GCGGCCGC",
    }
    seq = clean_seq(seq)
    return ";".join([name for name, site in sites.items() if site in seq])


def rank_pct(s):
    return pd.to_numeric(s, errors="coerce").rank(pct=True)


def optional_bool_gate(df, col, default=True):
    if col not in df.columns:
        return pd.Series(default, index=df.index)
    s = df[col]
    if pd.api.types.is_bool_dtype(s):
        return s.fillna(False)
    return s.astype(str).str.strip().str.lower().isin(["true", "1", "yes", "y"])


def nonempty_series(df, col):
    if col not in df.columns:
        return pd.Series(dtype=str)
    vals = df[col].dropna().astype(str).str.strip()
    return vals[(vals != "") & (vals.str.lower() != "nan")]


def choose_gene_column(df):
    for col in ["gene_name", "gene_id"]:
        if len(nonempty_series(df, col)):
            return col
    return None


def write_final_diversity_summary(lib, path):
    gene_col = choose_gene_column(lib)
    genes = nonempty_series(lib, gene_col) if gene_col else pd.Series(dtype=str)
    seq_clusters = nonempty_series(lib, "seq_cluster_id")
    gene_counts = genes.value_counts()
    seq_cluster_counts = seq_clusters.value_counts()

    lines = [
        "Final library gene and 5'UTR sequence-similarity cluster summary",
        "=" * 100,
        "Note: seq_cluster_id is a 5'UTR sequence-similarity cluster, not a gene cluster.",
        f"selected_n: {len(lib)}",
        f"gene_key_column: {gene_col or 'MISSING'}",
        f"n_unique_seq_clusters: {int(seq_cluster_counts.shape[0])}",
        f"max_per_seq_cluster: {int(seq_cluster_counts.max()) if len(seq_cluster_counts) else 0}",
        f"n_unique_genes: {int(gene_counts.shape[0])}",
        f"max_per_gene: {int(gene_counts.max()) if len(gene_counts) else 0}",
        "",
        "[Top genes by selected candidate count]",
        gene_counts.head(25).to_string() if len(gene_counts) else "No gene_name/gene_id values available.",
        "",
        "[Top 5'UTR sequence-similarity clusters by selected candidate count]",
        seq_cluster_counts.head(25).to_string() if len(seq_cluster_counts) else "No seq_cluster_id values available.",
    ]
    path.write_text("\n".join(lines), encoding="utf-8")


def metric_summary(df, base_candidate_pool=None, evidence_candidate_pool=None, unfilled_n=0):
    gene_col = choose_gene_column(df)
    genes = nonempty_series(df, gene_col) if gene_col else pd.Series(dtype=str)
    clusters = nonempty_series(df, "seq_cluster_id")
    heavy = pd.to_numeric(df["heavy_ensemble_score"], errors="coerce") if "heavy_ensemble_score" in df.columns else pd.Series(np.nan, index=df.index)
    robust = pd.to_numeric(df["robust_public_te_rank"], errors="coerce") if "robust_public_te_rank" in df.columns else pd.Series(np.nan, index=df.index)
    return {
        "selected_n": len(df),
        "unfilled_n": int(unfilled_n),
        "base_candidate_pool": int(base_candidate_pool) if base_candidate_pool is not None else np.nan,
        "evidence_candidate_pool": int(evidence_candidate_pool) if evidence_candidate_pool is not None else np.nan,
        "heavy_ensemble_score_non_null": int(heavy.notna().sum()),
        "mean_heavy_ensemble_score": float(heavy.mean()) if heavy.notna().any() else np.nan,
        "mean_robust_public_te_rank": float(robust.mean()) if robust.notna().any() else np.nan,
        "n_unique_seq_clusters": int(clusters.nunique()),
        "max_per_seq_cluster": int(clusters.value_counts().max()) if len(clusters) else 0,
        "n_unique_genes": int(genes.nunique()),
        "max_per_gene": int(genes.value_counts().max()) if len(genes) else 0,
    }


def write_uaug0_validation_report(lib, path, requested_n, max_cluster_cap, max_gene_cap):
    gene_col = choose_gene_column(lib)
    genes = nonempty_series(lib, gene_col) if gene_col else pd.Series(dtype=str)
    clusters = nonempty_series(lib, "seq_cluster_id")
    heavy = pd.to_numeric(lib["heavy_ensemble_score"], errors="coerce") if "heavy_ensemble_score" in lib.columns else pd.Series(np.nan, index=lib.index)
    robust = pd.to_numeric(lib["robust_public_te_rank"], errors="coerce") if "robust_public_te_rank" in lib.columns else pd.Series(np.nan, index=lib.index)
    uaug = pd.to_numeric(lib["uaug_count"], errors="coerce") if "uaug_count" in lib.columns else pd.Series(np.nan, index=lib.index)
    max_per_cluster = int(clusters.value_counts().max()) if len(clusters) else 0
    max_per_gene = int(genes.value_counts().max()) if len(genes) else 0
    uaug_positive = int(uaug.fillna(999).gt(0).sum()) if len(uaug) else 0
    lines = [
        "uAUG=0 production validation report",
        "=" * 100,
        f"selected_n: {len(lib)}",
        f"requested_n: {requested_n}",
        f"uaug_positive_n: {uaug_positive}",
        f"uaug0_policy_pass: {uaug_positive == 0}",
        "",
        "[Cluster diversity]",
        f"n_unique_seq_clusters: {int(clusters.nunique())}",
        f"max_per_seq_cluster: {max_per_cluster}",
        f"cluster_cap: {max_cluster_cap}",
        f"cluster_cap_pass: {max_per_cluster <= max_cluster_cap}",
        "",
        "[Gene diversity]",
        f"gene_column: {gene_col or 'NA'}",
        f"n_unique_genes: {int(genes.nunique())}",
        f"max_per_gene: {max_per_gene}",
        f"gene_cap: {max_gene_cap}",
        f"gene_cap_pass: {max_per_gene <= max_gene_cap}",
        "",
        "[Evidence means]",
        f"mean_heavy_ensemble_score: {float(heavy.mean()) if heavy.notna().any() else 'NA'}",
        f"mean_robust_public_te_rank: {float(robust.mean()) if robust.notna().any() else 'NA'}",
    ]
    path.write_text("\n".join(lines), encoding="utf-8")


def write_uaug_audit(lib):
    if "uaug_count" not in lib.columns:
        return
    x = lib.copy()
    x["uaug_count_numeric"] = pd.to_numeric(x["uaug_count"], errors="coerce")
    group_cols = [c for c in ["library_group", "selection_source"] if c in x.columns]
    rows = []
    for key, g in x.groupby(group_cols, dropna=False) if group_cols else [("all", x)]:
        if not isinstance(key, tuple):
            key = (key,)
        row = {col: val for col, val in zip(group_cols, key)}
        uaug_pos = g["uaug_count_numeric"].fillna(999).gt(0)
        heavy = pd.to_numeric(g["heavy_ensemble_score"], errors="coerce") if "heavy_ensemble_score" in g.columns else pd.Series(np.nan, index=g.index)
        robust = pd.to_numeric(g["robust_public_te_rank"], errors="coerce") if "robust_public_te_rank" in g.columns else pd.Series(np.nan, index=g.index)
        clusters = nonempty_series(g, "seq_cluster_id")
        row.update({
            "selected_n": len(g),
            "uaug0_n": int((g["uaug_count_numeric"] == 0).sum()),
            "uaug_positive_n": int(uaug_pos.sum()),
            "uaug_positive_ratio": float(uaug_pos.mean()) if len(g) else np.nan,
            "mean_heavy_ensemble_score": float(heavy.mean()) if heavy.notna().any() else np.nan,
            "mean_robust_public_te_rank": float(robust.mean()) if robust.notna().any() else np.nan,
            "n_unique_seq_clusters": int(clusters.nunique()),
            "max_per_seq_cluster": int(clusters.value_counts().max()) if len(clusters) else 0,
        })
        rows.append(row)
    summary = pd.DataFrame(rows)
    summary.to_csv(OUT_UAUG_SUMMARY_CSV, index=False)

    positive_cols = [c for c in [
        "library_index", "library_group", "selection_source", "utr_id", "gene_id", "gene_name",
        "seq_cluster_id", "uaug_count", "length", "gc_content", "robust_public_te_rank",
        "heavy_ensemble_score", "is_expressed_public", "expression_qc_reason", SEQ,
    ] if c in x.columns]
    x[x["uaug_count_numeric"].fillna(0).gt(0)][positive_cols].to_csv(OUT_UAUG_POSITIVE_CSV, index=False)

    total_pos = int(x["uaug_count_numeric"].fillna(0).gt(0).sum())
    lines = [
        "uAUG source audit for final selected library",
        "=" * 100,
        f"selected_n: {len(x)}",
        f"uaug0_n: {int((x['uaug_count_numeric'] == 0).sum())}",
        f"uaug_positive_n: {total_pos}",
        f"uaug_positive_ratio: {total_pos / len(x):.4f}" if len(x) else "uaug_positive_ratio: NA",
        "",
        "[By library_group and selection_source]",
        summary.to_string(index=False) if len(summary) else "No summary rows.",
        "",
        f"Saved CSV summary: {OUT_UAUG_SUMMARY_CSV}",
        f"Saved positive rows: {OUT_UAUG_POSITIVE_CSV}",
    ]
    OUT_UAUG_SUMMARY_TXT.write_text("\n".join(lines), encoding="utf-8")


def dry_run_select_uaug0(base_cand, args, quotas, production_lib):
    uaug0 = pd.to_numeric(base_cand["uaug_count"], errors="coerce").fillna(999).eq(0)
    base0 = base_cand[uaug0].copy().reset_index(drop=True)
    evidence0 = base0[
        optional_bool_gate(base0, "is_expressed_public", default=True) &
        pd.to_numeric(base0["robust_public_te_rank"], errors="coerce").notna()
    ].copy()

    selected = []
    used_seq = set()
    cluster_counts = defaultdict(int)
    taken_by_group = defaultdict(int)

    def can_take(row, max_per):
        seq = row[SEQ]
        cid = str(row.get("seq_cluster_id", seq))
        return seq not in used_seq and cluster_counts[cid] < max_per

    def take(pool, n, group, sort_cols, ascending=None, max_per_cluster=None, source=None):
        source = source or group
        max_per_cluster = args.max_per_cluster if max_per_cluster is None else max_per_cluster
        ascending = [False] * len(sort_cols) if ascending is None else ascending
        rows = []
        p = pool.copy().sort_values(sort_cols, ascending=ascending)
        for _, row in p.iterrows():
            if len(rows) >= n:
                break
            if can_take(row, max_per_cluster):
                rows.append(row)
                used_seq.add(row[SEQ])
                cluster_counts[str(row.get("seq_cluster_id", row[SEQ]))] += 1
        if rows:
            out = pd.DataFrame(rows)
            out["library_group"] = group
            out["selection_source"] = source
            selected.append(out)
        taken_by_group[group] += len(rows)
        return len(rows)

    take(evidence0, quotas["A_publicTE_high_confidence"], "A_publicTE_high_confidence", ["robust_public_te_rank", "day_consensus_TE_rank" if "day_consensus_TE_rank" in evidence0.columns else "cluster_diverse_evidence_score"], source="uaug0_evidence_cand")
    take(evidence0, quotas["B_TE_model_classifier_supported"], "B_TE_model_classifier_supported", ["model_support_score", "robust_public_te_rank"], source="uaug0_evidence_cand")
    if "protein_abundance_rank" in evidence0.columns:
        take(evidence0[evidence0["protein_abundance_rank"].notna()], quotas["C_protein_abundance_supported"], "C_protein_abundance_supported", ["protein_abundance_rank", "cluster_diverse_evidence_score"], source="uaug0_evidence_cand")
    if "protein_residual_rank" in evidence0.columns:
        take(evidence0[evidence0["protein_residual_rank"].notna()], quotas["D_protein_residual_supported"], "D_protein_residual_supported", ["protein_residual_rank", "cluster_diverse_evidence_score"], source="uaug0_evidence_cand")
    if "multi_omics_utr_rank" in evidence0.columns:
        take(evidence0, quotas["E_multiomics_consensus_high"], "E_multiomics_consensus_high", ["multi_omics_utr_rank", "cluster_diverse_evidence_score"], source="uaug0_evidence_cand")

    robust_rank = pd.to_numeric(base0["robust_public_te_rank"], errors="coerce")
    exploratory = base0[robust_rank.between(0.35, 0.90) | robust_rank.isna()].copy()
    if len(exploratory):
        exploratory["length_bin"] = pd.cut(exploratory["length"], [49, 60, 75, 90, 100], labels=["50-60", "61-75", "76-90", "91-100"])
        exploratory["gc_bin"] = pd.cut(exploratory["gc_content"], [0.299, 0.45, 0.60, 0.75], labels=["30-45", "45-60", "60-75"])
        exploratory["diversity_key"] = exploratory["length_bin"].astype(str) + ":" + exploratory["gc_bin"].astype(str) + ":uAUG" + exploratory["uaug_count"].astype(str)
        exp_rows = []
        for _, g in exploratory.sort_values("cluster_diverse_evidence_score", ascending=False).groupby("diversity_key", dropna=False):
            exp_rows.append(g.head(max(1, quotas["F_sequence_diverse_exploratory"] // max(1, exploratory["diversity_key"].nunique()))))
        exp_pool = pd.concat(exp_rows).sort_values("cluster_diverse_evidence_score", ascending=False) if exp_rows else exploratory
    else:
        exp_pool = exploratory
    take(exp_pool, quotas["F_sequence_diverse_exploratory"], "F_sequence_diverse_exploratory", ["cluster_diverse_evidence_score", "robust_public_te_rank"], source="uaug0_base_cand_exploratory")

    diversity = base0.copy()
    diversity["diversity_score"] = (
        (1 - abs(pd.to_numeric(diversity["gc_content"], errors="coerce") - 0.52)) +
        (1 - abs(pd.to_numeric(diversity["length"], errors="coerce") - 75) / 50)
    )
    take(diversity, quotas["G_length_GC_uAUG_diversity"], "G_length_GC_uAUG_diversity", ["diversity_score", "cluster_diverse_evidence_score"], source="uaug0_base_cand_diversity")

    low = base0.copy()
    low["low_score"] = 1 - pd.to_numeric(low["robust_public_te_rank"], errors="coerce").fillna(1)
    take(low[pd.to_numeric(low["robust_public_te_rank"], errors="coerce") <= 0.30], quotas["H_low_signal_negative_controls"], "H_low_signal_negative_controls", ["low_score", "gc_content"], source="uaug0_base_cand_low_publicTE")
    if taken_by_group["H_low_signal_negative_controls"] < quotas["H_low_signal_negative_controls"]:
        no_evidence_low = low[pd.to_numeric(low["robust_public_te_rank"], errors="coerce").isna()].copy()
        no_evidence_low["clean_control_score"] = (
            (1 - abs(pd.to_numeric(no_evidence_low["gc_content"], errors="coerce") - 0.52)) +
            (1 - abs(pd.to_numeric(no_evidence_low["length"], errors="coerce") - 75) / 50)
        )
        take(no_evidence_low, quotas["H_low_signal_negative_controls"] - taken_by_group["H_low_signal_negative_controls"], "H_low_signal_negative_controls", ["clean_control_score"], source="uaug0_base_cand_no_TE_clean_control")

    lib = pd.concat(selected, ignore_index=True) if selected else base0.head(0)
    if len(lib) < args.n:
        used_seq = set(lib[SEQ])
        cluster_counts = defaultdict(int, lib.groupby("seq_cluster_id").size().to_dict() if "seq_cluster_id" in lib.columns else {})

        def fill_from(pool, source):
            nonlocal lib
            rows = []
            p = pool[~pool[SEQ].isin(used_seq)].sort_values(["cluster_diverse_evidence_score", "robust_public_te_rank"], ascending=False)
            for _, row in p.iterrows():
                if len(lib) + len(rows) >= args.n:
                    break
                cid = str(row.get("seq_cluster_id", row[SEQ]))
                if cluster_counts[cid] < args.allow_cluster_fill:
                    row = row.copy()
                    row["library_group"] = "A_publicTE_high_confidence"
                    row["selection_source"] = source
                    row["selection_phase"] = "evidence_refill"
                    rows.append(row)
                    used_seq.add(row[SEQ])
                    cluster_counts[cid] += 1
            if rows:
                lib = pd.concat([lib, pd.DataFrame(rows)], ignore_index=True, sort=False)
                taken_by_group["A_publicTE_high_confidence"] += len(rows)

        fill_from(evidence0, "uaug0_fill_evidence_cand")

    lib = lib.drop_duplicates(subset=[SEQ], keep="first").head(args.n).copy()
    lib["library_index"] = np.arange(1, len(lib) + 1)
    unfilled_n = max(0, args.n - len(lib))

    summary = metric_summary(lib, base_candidate_pool=len(base0), evidence_candidate_pool=len(evidence0), unfilled_n=unfilled_n)
    summary.update({
        "mode": "uaug0_hard_filter_dry_run",
        "selected_uaug_positive_n": int(pd.to_numeric(lib.get("uaug_count", pd.Series(dtype=float)), errors="coerce").fillna(0).gt(0).sum()),
        "selected_uaug0_n": int(pd.to_numeric(lib.get("uaug_count", pd.Series(dtype=float)), errors="coerce").fillna(999).eq(0).sum()),
    })
    pd.DataFrame([summary]).to_csv(OUT_UAUG0_DRY_SUMMARY, index=False)

    shortfall_rows = []
    for group, quota in quotas.items():
        selected_n = int(taken_by_group[group])
        shortfall_rows.append({
            "library_group": group,
            "quota": quota,
            "selected_uaug0_dry_run": selected_n,
            "shortfall": max(0, quota - selected_n),
            "source_pool": "uaug0_evidence_cand" if group.startswith(("A_", "B_", "C_", "D_", "E_")) else "uaug0_base_cand",
            "reason": "insufficient uaug0 candidates under cluster/sequence constraints" if selected_n < quota else "filled",
        })
    pd.DataFrame(shortfall_rows).to_csv(OUT_UAUG0_SHORTFALL, index=False)

    prod_seq = set(production_lib[SEQ]) if SEQ in production_lib.columns else set()
    replacement_cols = [c for c in [
        "library_index", "library_group", "selection_source", "utr_id", "gene_id", "gene_name",
        "seq_cluster_id", "uaug_count", "length", "gc_content", "robust_public_te_rank",
        "heavy_ensemble_score", "is_expressed_public", "expression_qc_reason", SEQ,
    ] if c in lib.columns]
    lib[~lib[SEQ].isin(prod_seq)][replacement_cols].to_csv(OUT_UAUG0_REPLACEMENTS, index=False)


def prep(df):
    df = df.copy()
    df[SEQ] = df[SEQ].map(clean_seq)
    if "utr5_length_final" in df.columns:
        df["length"] = pd.to_numeric(df["utr5_length_final"], errors="coerce")
    elif "utr5_length_tss_corrected" in df.columns:
        df["length"] = pd.to_numeric(df["utr5_length_tss_corrected"], errors="coerce")
    else:
        df["length"] = df[SEQ].str.len()
    if "gc_content" not in df.columns:
        df["gc_content"] = df[SEQ].map(gc_content)
    if "uaug_count" not in df.columns:
        df["uaug_count"] = df[SEQ].str.count("ATG")
    if "seq_cluster_id" not in df.columns:
        df["seq_cluster_id"] = "NOCLUSTER_" + df.index.astype(str)
    if "has_proteomics_label" not in df.columns:
        prot_cols = [c for c in ["protein_abundance_rank", "protein_residual_rank"] if c in df.columns]
        df["has_proteomics_label"] = df[prot_cols].notna().any(axis=1) if prot_cols else False
    df["forbidden_sites"] = df[SEQ].map(forbidden_sites)
    return df


def evidence_score(df):
    z = pd.Series(0.0, index=df.index)
    def val(c, fallback=0.0):
        if c in df.columns:
            return pd.to_numeric(df[c], errors="coerce").fillna(fallback)
        return pd.Series(fallback, index=df.index)

    model_score = pd.Series(0.0, index=df.index)
    model_parts = []
    for c in ["heavy_ensemble_score", "automl_ensemble_score", "proteomics_enriched_score", "model_pred_rank_40_200train"]:
        if c in df.columns:
            model_parts.append(val(c, 0.0))
    if model_parts:
        model_score = sum(model_parts) / len(model_parts)

    z = (
        0.42 * val("robust_public_te_rank") +
        0.16 * val("day_consensus_TE_rank") +
        0.12 * val("protein_residual_rank", val("robust_public_te_rank")) +
        0.10 * val("protein_abundance_rank", val("robust_public_te_rank")) +
        0.10 * val("multi_omics_utr_rank", val("robust_public_te_rank")) +
        0.07 * model_score +
        0.03 * val("tss_confidence_score")
    )
    return z


def parse_reference_controls():
    out = []
    csv_path = BASE / "00_raw_data/04_reference_controls/reference_controls.csv"
    if csv_path.exists():
        try:
            rc = pd.read_csv(csv_path)
            seq_col = next((c for c in rc.columns if c.lower() in ["sequence", "seq", "utr5_sequence"]), None)
            if seq_col:
                for _, r in rc.iterrows():
                    seq = clean_seq(r[seq_col])
                    if seq:
                        out.append({
                            "control_id": str(r.get("control_id", r.get("id", f"control_{len(out)+1}"))),
                            "control_type": str(r.get("control_type", "reference_control")),
                            SEQ: seq,
                            "note": str(r.get("note", "")),
                        })
        except Exception as e:
            print("[WARN] failed reference_controls.csv", e)
    fasta_dir = BASE / "00_raw_data/04_reference_controls"
    for fa in list(fasta_dir.glob("*.fa")) + list(fasta_dir.glob("*.fasta")):
        cid = None
        seqs = []
        for line in fa.read_text(errors="ignore").splitlines():
            if line.startswith(">"):
                if cid and seqs:
                    out.append({"control_id": cid, "control_type": "fasta_control", SEQ: clean_seq("".join(seqs)), "note": fa.name})
                cid = line[1:].strip().split()[0]
                seqs = []
            else:
                seqs.append(line.strip())
        if cid and seqs:
            out.append({"control_id": cid, "control_type": "fasta_control", SEQ: clean_seq("".join(seqs)), "note": fa.name})
    return pd.DataFrame(out)


def write_fasta(df, path):
    with open(path, "w", encoding="utf-8") as out:
        for _, r in df.iterrows():
            idx = r.get("library_index", "NA")
            gene = str(r.get("gene_name", r.get("control_id", "NA"))).replace(" ", "_").replace("/", "_")
            group = str(r.get("library_group", "NA")).replace(" ", "_")
            seq = clean_seq(r[SEQ])
            out.write(f">lib{idx}|{group}|{gene}|len={len(seq)}\n")
            for i in range(0, len(seq), 80):
                out.write(seq[i:i+80] + "\n")


def main():
    ap = argparse.ArgumentParser(description="Select 2000 cluster-diverse, evidence-balanced 50-100bp CHO 5'UTR library")
    ap.add_argument("--input", default=None)
    ap.add_argument("--n", type=int, default=2000)
    ap.add_argument("--max-per-cluster", type=int, default=1)
    ap.add_argument("--allow-cluster-fill", type=int, default=2, help="allow this per cluster during final fill if needed")
    ap.add_argument("--max-per-gene", type=int, default=3, help="maximum selected candidates per gene_name/gene_id when available")
    ap.add_argument(
        "--allow-evidence-shortfall",
        action="store_true",
        help="write a library smaller than --n instead of failing when evidence-only refill is insufficient",
    )
    args = ap.parse_args()

    path = choose_input(args.input)
    print("[LOAD]", path)
    df = prep(pd.read_csv(path))
    expression_gate = optional_bool_gate(df, "is_expressed_public", default=True)
    base_cand = df[
        df["length"].between(50, 100) &
        pd.to_numeric(df["gc_content"], errors="coerce").between(0.30, 0.75) &
        (pd.to_numeric(df["uaug_count"], errors="coerce").fillna(999) == 0) &
        df["forbidden_sites"].eq("") &
        df[SEQ].str.len().gt(0) &
        ~df[SEQ].str.contains("N", regex=False)
    ].copy()
    base_cand = base_cand.drop_duplicates(subset=[SEQ]).reset_index(drop=True)
    if "robust_public_te_rank" not in base_cand.columns:
        base_cand["robust_public_te_rank"] = np.nan
    base_cand["cluster_diverse_evidence_score"] = evidence_score(base_cand)
    base_expression_gate = optional_bool_gate(base_cand, "is_expressed_public", default=True)
    evidence_cand = base_cand[
        base_expression_gate &
        pd.to_numeric(base_cand["robust_public_te_rank"], errors="coerce").notna()
    ].copy()

    print(f"[BASE CANDIDATES] {len(base_cand):,}")
    print(f"[EVIDENCE CANDIDATES] {len(evidence_cand):,}")
    if len(base_cand) < args.n:
        print(f"[WARN] base candidate count {len(base_cand)} < requested {args.n}; output will be smaller unless controls add rows")

    selected = []
    used_seq = set()
    cluster_counts = defaultdict(int)
    gene_counts = defaultdict(int)
    gene_col_for_cap = choose_gene_column(base_cand)

    def gene_key(row):
        if not gene_col_for_cap:
            return None
        val = row.get(gene_col_for_cap)
        if pd.isna(val) or not str(val).strip():
            return None
        return str(val)

    def can_take(row, max_per):
        seq = row[SEQ]
        cid = str(row.get("seq_cluster_id", seq))
        gkey = gene_key(row)
        gene_ok = gkey is None or gene_counts[gkey] < args.max_per_gene
        return seq not in used_seq and cluster_counts[cid] < max_per and gene_ok

    def take(pool, n, group, sort_cols, ascending=None, max_per_cluster=None, source=None):
        nonlocal selected, used_seq, cluster_counts
        source = source or group
        if max_per_cluster is None:
            max_per_cluster = args.max_per_cluster
        if ascending is None:
            ascending = [False] * len(sort_cols)
        p = pool.copy()
        p = p.sort_values(sort_cols, ascending=ascending)
        rows = []
        for _, row in p.iterrows():
            if len(rows) >= n:
                break
            if can_take(row, max_per_cluster):
                rows.append(row)
                used_seq.add(row[SEQ])
                cluster_counts[str(row.get("seq_cluster_id", row[SEQ]))] += 1
                gkey = gene_key(row)
                if gkey is not None:
                    gene_counts[gkey] += 1
        if rows:
            out = pd.DataFrame(rows)
            out["library_group"] = group
            out["selection_source"] = source
            selected.append(out)
            print(f"  take {group}: {len(out)}/{n}")
        else:
            print(f"  take {group}: 0/{n}")
        return len(rows)

    # v1.4 primary targets deliberately differ from v1.3: expand non-J evidence
    # groups and reduce H controls. Any remaining shortage is filled by K1-K4.
    quotas = {
        "A_publicTE_high_confidence": 550,
        "B_TE_model_classifier_supported": 350,
        "C_protein_abundance_supported": 300,
        "D_protein_residual_supported": 250,
        "E_multiomics_consensus_high": 300,
        "F_sequence_diverse_exploratory": 150,
        "G_length_GC_uAUG_diversity": 50,
        "H_low_signal_negative_controls": 50,
    }

    model_cols = [c for c in ["heavy_ensemble_score", "automl_ensemble_score", "proteomics_enriched_score", "model_pred_rank_40_200train"] if c in base_cand.columns]
    if model_cols:
        numeric_model = base_cand[model_cols].apply(pd.to_numeric, errors="coerce")
        base_cand["model_support_score"] = numeric_model.mean(axis=1)
        base_cand["has_classifier_support"] = numeric_model.notna().any(axis=1)
    else:
        base_cand["model_support_score"] = base_cand["cluster_diverse_evidence_score"]
        base_cand["has_classifier_support"] = False
    protein_cols = [c for c in ["protein_abundance_rank", "protein_residual_rank"] if c in base_cand.columns]
    protein_numeric = (
        base_cand[protein_cols].apply(pd.to_numeric, errors="coerce").notna().any(axis=1)
        if protein_cols else pd.Series(False, index=base_cand.index)
    )
    base_cand["has_protein_support"] = (
        optional_bool_gate(base_cand, "has_proteomics_label", default=False) | protein_numeric
    )
    base_cand["has_multiomics_support"] = (
        pd.to_numeric(base_cand["multi_omics_utr_rank"], errors="coerce").notna()
        if "multi_omics_utr_rank" in base_cand.columns
        else False
    )
    evidence_cand = base_cand[
        base_expression_gate &
        pd.to_numeric(base_cand["robust_public_te_rank"], errors="coerce").notna()
    ].copy()
    take(evidence_cand, quotas["A_publicTE_high_confidence"], "A_publicTE_high_confidence", ["robust_public_te_rank", "day_consensus_TE_rank" if "day_consensus_TE_rank" in evidence_cand.columns else "cluster_diverse_evidence_score"], source="evidence_cand")
    take(
        evidence_cand[evidence_cand["has_classifier_support"]],
        quotas["B_TE_model_classifier_supported"],
        "B_TE_model_classifier_supported",
        ["model_support_score", "robust_public_te_rank"],
        source="evidence_cand_classifier_supported",
    )

    if "protein_abundance_rank" in evidence_cand.columns:
        take(evidence_cand[evidence_cand["protein_abundance_rank"].notna()], quotas["C_protein_abundance_supported"], "C_protein_abundance_supported", ["protein_abundance_rank", "cluster_diverse_evidence_score"], source="evidence_cand")
    if "protein_residual_rank" in evidence_cand.columns:
        take(evidence_cand[evidence_cand["protein_residual_rank"].notna()], quotas["D_protein_residual_supported"], "D_protein_residual_supported", ["protein_residual_rank", "cluster_diverse_evidence_score"], source="evidence_cand")
    if "multi_omics_utr_rank" in evidence_cand.columns:
        take(
            evidence_cand[evidence_cand["has_multiomics_support"]],
            quotas["E_multiomics_consensus_high"],
            "E_multiomics_consensus_high",
            ["multi_omics_utr_rank", "cluster_diverse_evidence_score"],
            source="evidence_cand_multiomics_supported",
        )

    # Exploratory: start with mid/high public TE when available, then allow clean no-evidence sequences.
    robust_rank = pd.to_numeric(base_cand["robust_public_te_rank"], errors="coerce")
    exploratory = base_cand[robust_rank.between(0.35, 0.90)].copy()
    no_evidence_exploratory = base_cand[robust_rank.isna()].copy()
    if len(no_evidence_exploratory):
        exploratory = pd.concat([exploratory, no_evidence_exploratory], ignore_index=True, sort=False)
    exploratory["length_bin"] = pd.cut(exploratory["length"], [49, 60, 75, 90, 100], labels=["50-60", "61-75", "76-90", "91-100"])
    exploratory["gc_bin"] = pd.cut(exploratory["gc_content"], [0.299, 0.45, 0.60, 0.75], labels=["30-45", "45-60", "60-75"])
    exploratory["diversity_key"] = exploratory["length_bin"].astype(str) + ":" + exploratory["gc_bin"].astype(str) + ":uAUG" + exploratory["uaug_count"].astype(str)
    # Interleave by diversity key
    exp_rows = []
    for _, g in exploratory.sort_values("cluster_diverse_evidence_score", ascending=False).groupby("diversity_key", dropna=False):
        exp_rows.append(g.head(max(1, quotas["F_sequence_diverse_exploratory"] // max(1, exploratory["diversity_key"].nunique()))))
    exp_pool = pd.concat(exp_rows).sort_values("cluster_diverse_evidence_score", ascending=False) if exp_rows else exploratory
    take(exp_pool, quotas["F_sequence_diverse_exploratory"], "F_sequence_diverse_exploratory", ["cluster_diverse_evidence_score", "robust_public_te_rank"], source="base_cand_exploratory")

    diversity = base_cand.copy()
    diversity["diversity_score"] = (
        (1 - abs(pd.to_numeric(diversity["gc_content"], errors="coerce") - 0.52)) +
        (1 - abs(pd.to_numeric(diversity["length"], errors="coerce") - 75) / 50) +
        0.2 * (diversity["uaug_count"] == 0).astype(float)
    )
    take(diversity, quotas["G_length_GC_uAUG_diversity"], "G_length_GC_uAUG_diversity", ["diversity_score", "cluster_diverse_evidence_score"], source="base_cand_diversity")

    low = base_cand.copy()
    low["low_score"] = 1 - pd.to_numeric(low["robust_public_te_rank"], errors="coerce").fillna(1)
    if "protein_abundance_rank" in low.columns:
        low["low_score"] += (1 - pd.to_numeric(low["protein_abundance_rank"], errors="coerce").fillna(1)) * 0.25
    if "protein_residual_rank" in low.columns:
        low["low_score"] += (1 - pd.to_numeric(low["protein_residual_rank"], errors="coerce").fillna(1)) * 0.25
    h_taken = take(low[pd.to_numeric(low["robust_public_te_rank"], errors="coerce") <= 0.30], quotas["H_low_signal_negative_controls"], "H_low_signal_negative_controls", ["low_score", "gc_content"], source="base_cand_low_publicTE")
    if h_taken < quotas["H_low_signal_negative_controls"]:
        no_evidence_low = low[pd.to_numeric(low["robust_public_te_rank"], errors="coerce").isna()].copy()
        no_evidence_low["clean_control_score"] = (
            (1 - abs(pd.to_numeric(no_evidence_low["gc_content"], errors="coerce") - 0.52)) +
            (1 - abs(pd.to_numeric(no_evidence_low["length"], errors="coerce") - 75) / 50) +
            0.2 * (pd.to_numeric(no_evidence_low["uaug_count"], errors="coerce").fillna(999) == 0).astype(float)
        )
        take(no_evidence_low, quotas["H_low_signal_negative_controls"] - h_taken, "H_low_signal_negative_controls", ["clean_control_score"], source="base_cand_no_TE_clean_control")

    lib = pd.concat(selected, ignore_index=True) if selected else base_cand.head(0)

    # add reference controls up to 50 if present, without enforcing cluster.
    controls = parse_reference_controls()
    if len(controls):
        controls = prep(controls)
        controls = controls[pd.to_numeric(controls["uaug_count"], errors="coerce").fillna(999).eq(0)].copy()
        controls["library_group"] = "I_reference_controls"
        controls["selection_source"] = "reference_control"
        controls["is_reference_control"] = True
        controls["cluster_diverse_evidence_score"] = np.nan
        controls = controls.head(50)
        lib["is_reference_control"] = False
        lib = pd.concat([lib, controls], ignore_index=True, sort=False)
        print(f"  added reference controls: {len(controls)}")
    else:
        lib["is_reference_control"] = False
        print("  reference controls: none found")

    # Ordered non-J refill. Each K pool is constrained by the same sequence QC,
    # gene cap, and relaxed sequence-cluster cap used for the final library.
    refill_audit_rows = []
    if len(lib) < args.n:
        used_seq = set(lib[SEQ])
        cluster_counts = defaultdict(int, lib[~lib["is_reference_control"].fillna(False)].groupby("seq_cluster_id").size().to_dict() if "seq_cluster_id" in lib.columns else {})
        gene_col_for_cap = choose_gene_column(lib)
        gene_counts = defaultdict(int, nonempty_series(lib, gene_col_for_cap).value_counts().to_dict() if gene_col_for_cap else {})
        fill = []
        fill_group_counts = defaultdict(int)

        def fill_from(pool, group, source, sort_cols, limit=None):
            nonlocal fill
            available_sort_cols = [c for c in sort_cols if c in pool.columns]
            if not available_sort_cols:
                available_sort_cols = ["cluster_diverse_evidence_score"]
            p = pool[~pool[SEQ].isin(used_seq)].sort_values(available_sort_cols, ascending=False)
            selected_this_call = 0
            for _, row in p.iterrows():
                if len(lib) + len(fill) >= args.n:
                    break
                if limit is not None and selected_this_call >= limit:
                    break
                cid = str(row.get("seq_cluster_id", row[SEQ]))
                gkey = gene_key(row)
                gene_ok = gkey is None or gene_counts[gkey] < args.max_per_gene
                if cluster_counts[cid] < args.allow_cluster_fill and gene_ok:
                    row = row.copy()
                    row["selection_source"] = source
                    row["selection_phase"] = "non_j_relaxed_refill"
                    row["library_group"] = group
                    fill.append(row)
                    used_seq.add(row[SEQ])
                    cluster_counts[cid] += 1
                    if gkey is not None:
                        gene_counts[gkey] += 1
                    fill_group_counts[group] += 1
                    selected_this_call += 1
            return selected_this_call

        k4_pool = base_cand.copy()
        k4_pool["length_bin"] = pd.cut(
            pd.to_numeric(k4_pool["length"], errors="coerce"),
            [49, 60, 75, 90, 100],
            labels=False,
        )
        k4_pool["gc_bin"] = pd.cut(
            pd.to_numeric(k4_pool["gc_content"], errors="coerce"),
            [0.299, 0.45, 0.60, 0.75],
            labels=False,
        )
        diversity_frequency = k4_pool.groupby(["length_bin", "gc_bin"], dropna=False)[SEQ].transform("count")
        k4_pool["relaxed_diversity_score"] = (
            1.0 / diversity_frequency.clip(lower=1) +
            0.25 * pd.to_numeric(k4_pool["cluster_diverse_evidence_score"], errors="coerce").fillna(0)
        )
        refill_specs = [
            (
                "K1_ABE_evidence_relaxed",
                evidence_cand,
                "K1_unselected_ABE_evidence",
                ["robust_public_te_rank", "multi_omics_utr_rank", "cluster_diverse_evidence_score"],
            ),
            (
                "K2_CD_proteomics_relaxed",
                base_cand[base_cand["has_protein_support"]],
                "K2_unselected_CD_proteomics",
                ["protein_abundance_rank", "protein_residual_rank", "cluster_diverse_evidence_score"],
            ),
            (
                "K3_classifier_model_relaxed",
                base_cand[base_cand["has_classifier_support"]],
                "K3_unselected_classifier_model",
                ["model_support_score", "robust_public_te_rank", "cluster_diverse_evidence_score"],
            ),
            (
                "K4_FG_diversity_relaxed",
                k4_pool,
                "K4_unselected_FG_diversity",
                ["relaxed_diversity_score", "cluster_diverse_evidence_score"],
            ),
        ]
        for group, pool, source, sort_cols in refill_specs:
            target = args.n - len(lib) - len(fill)
            selected_n = fill_from(pool, group, source, sort_cols)
            refill_audit_rows.append({
                "library_group": group,
                "selection_source": source,
                "refill_pass": "ordered_K",
                "target_refill_n": target,
                "available_pool_n": len(pool),
                "selected_refill_n": selected_n,
            })
            if len(lib) + len(fill) >= args.n:
                break

        if fill:
            f = pd.DataFrame(fill)
            f["is_reference_control"] = False
            lib = pd.concat([lib, f], ignore_index=True, sort=False)
            print(f"  non-J K refill: {len(f)}; groups={dict(fill_group_counts)}")

    pd.DataFrame(
        refill_audit_rows,
        columns=[
            "library_group", "selection_source", "refill_pass", "target_refill_n",
            "available_pool_n", "selected_refill_n",
        ],
    ).to_csv(OUT_REFILL_AUDIT, index=False)

    lib = lib.drop_duplicates(subset=[SEQ], keep="first").head(args.n).copy()

    if len(lib) < args.n and not args.allow_evidence_shortfall:
        raise SystemExit(
            f"Non-J refill shortfall: selected {len(lib)} of {args.n}. "
            "J_fill is disabled and the filtered K1-K4 pools could not satisfy the caps. Review the "
            "candidate pool or use "
            "--allow-evidence-shortfall for audit-only output."
        )

    group_text = lib["library_group"].fillna("").astype(str)
    j_fill_selected_n = int(group_text.str.contains("J_fill", regex=False).sum())
    if j_fill_selected_n:
        raise SystemExit(f"J_fill policy violation: {j_fill_selected_n} J_fill rows selected")

    final_gene_col = choose_gene_column(lib)
    final_genes = nonempty_series(lib, final_gene_col) if final_gene_col else pd.Series(dtype=str)
    final_clusters = nonempty_series(lib, "seq_cluster_id")
    final_max_gene = int(final_genes.value_counts().max()) if len(final_genes) else 0
    final_max_cluster = int(final_clusters.value_counts().max()) if len(final_clusters) else 0
    if final_max_gene > args.max_per_gene:
        raise SystemExit(f"Gene cap violation: observed {final_max_gene}, allowed {args.max_per_gene}")
    if final_max_cluster > args.allow_cluster_fill:
        raise SystemExit(
            f"Sequence-cluster cap violation: observed {final_max_cluster}, allowed {args.allow_cluster_fill}"
        )

    prefix_counts = {
        prefix: int(group_text.str.startswith("K" if prefix == "K" else f"{prefix}_").sum())
        for prefix in "ABCDEFGHK"
    }
    protein_selected = int(optional_bool_gate(lib, "has_protein_support", default=False).sum())
    classifier_selected = int(optional_bool_gate(lib, "has_classifier_support", default=False).sum())
    multiomics_selected = int(optional_bool_gate(lib, "has_multiomics_support", default=False).sum())
    selection_qc = {
        "selected_n": len(lib),
        "requested_n": args.n,
        "shortage_n": max(0, args.n - len(lib)),
        "J_fill_selected_n": j_fill_selected_n,
        **{f"{prefix}_count": prefix_counts[prefix] for prefix in "ABCDEFGHK"},
        "total_protein_supported_selected_count": protein_selected,
        "total_classifier_supported_selected_count": classifier_selected,
        "total_multiomics_supported_selected_count": multiomics_selected,
        "max_per_gene": final_max_gene,
        "gene_cap": args.max_per_gene,
        "max_per_seq_cluster": final_max_cluster,
        "seq_cluster_cap": args.allow_cluster_fill,
    }
    pd.DataFrame([selection_qc]).to_csv(OUT_SELECTION_QC, index=False)

    lib["library_index"] = np.arange(1, len(lib) + 1)

    front = [c for c in [
        "library_index", "library_group", "is_reference_control", "utr_id", "gene_id", "gene_name", SEQ,
        "length", "gc_content", "uaug_count", "seq_cluster_id", "seq_cluster_size", "forbidden_sites",
        "selection_source", "is_expressed_public", "expression_qc_reason",
        "selection_phase",
        "robust_public_te_rank", "day_consensus_TE_rank", "protein_abundance_rank", "protein_residual_rank", "multi_omics_utr_rank",
        "cluster_diverse_evidence_score", "model_support_score", "has_proteomics_label",
        "has_protein_support", "has_classifier_support", "has_multiomics_support"
    ] if c in lib.columns]
    lib = lib[front + [c for c in lib.columns if c not in front]]
    lib.to_csv(OUT_CSV, index=False)
    write_fasta(lib, OUT_FASTA)
    write_final_diversity_summary(lib, OUT_DIVERSITY)
    write_uaug_audit(lib)
    write_uaug0_validation_report(lib, OUT_UAUG0_VALIDATION, args.n, args.allow_cluster_fill, args.max_per_gene)
    dry_run_select_uaug0(base_cand, args, quotas, lib)

    q = [
        "Cluster-diverse evidence-balanced 2000 library summary",
        "=" * 100,
        f"input: {path}",
        f"candidate_pool_after_QC: {len(base_cand)}",
        f"evidence_candidate_pool_after_expression_TE_QC: {len(evidence_cand)}",
        f"requested_n: {args.n}",
        f"selected_n: {len(lib)}",
        f"shortage_n: {selection_qc['shortage_n']}",
        f"max_per_cluster_primary: {args.max_per_cluster}",
        f"allow_cluster_fill: {args.allow_cluster_fill}",
        f"max_per_gene: {args.max_per_gene}",
        "uaug_policy: production_hard_filter_uaug_count_eq_0",
        "",
        "[Group counts]",
        lib["library_group"].value_counts(dropna=False).to_string(),
        "",
        "[Selection phase counts]",
        lib.get("selection_phase", pd.Series("primary", index=lib.index)).fillna("primary").value_counts().to_string(),
        "",
        "J_fill_policy: disabled",
        f"J_fill_selected_n: {selection_qc['J_fill_selected_n']}",
        "",
        "[Required v1.4 group totals]",
        "\n".join(f"{prefix}_count: {selection_qc[f'{prefix}_count']}" for prefix in "ABCDEFGHK"),
        f"total_protein_supported_selected_count: {selection_qc['total_protein_supported_selected_count']}",
        f"total_classifier_supported_selected_count: {selection_qc['total_classifier_supported_selected_count']}",
        f"total_multiomics_supported_selected_count: {selection_qc['total_multiomics_supported_selected_count']}",
        "",
        "[Selection source counts]",
        lib["selection_source"].value_counts(dropna=False).to_string() if "selection_source" in lib.columns else "NA",
        "",
        "[Cluster counts in final library]",
        lib["seq_cluster_id"].value_counts().describe().to_string() if "seq_cluster_id" in lib.columns else "No seq_cluster_id column",
        "",
        "[Length]",
        pd.to_numeric(lib["length"], errors="coerce").describe().to_string() if "length" in lib.columns else "NA",
        "",
        "[GC]",
        pd.to_numeric(lib["gc_content"], errors="coerce").describe().to_string() if "gc_content" in lib.columns else "NA",
        "",
        "[uAUG]",
        lib["uaug_count"].value_counts(dropna=False).sort_index().to_string() if "uaug_count" in lib.columns else "NA",
        "",
        "[Proteomics label coverage]",
        lib["has_proteomics_label"].value_counts(dropna=False).to_string() if "has_proteomics_label" in lib.columns else "NA",
        "",
        "[Public expression gate coverage]",
        (
            f"input_expression_pass: {int(expression_gate.sum())} / {len(df)} ({expression_gate.mean():.3f})\n"
            f"base_candidate_expression_pass: {int(base_expression_gate.sum())} / {len(base_cand)} ({base_expression_gate.mean():.3f})\n"
            f"evidence_candidates: {len(evidence_cand)}\n"
            + lib["is_expressed_public"].value_counts(dropna=False).to_string()
            if "is_expressed_public" in df.columns and "is_expressed_public" in lib.columns
            else "is_expressed_public column missing; expression gate not applied"
        ),
        "",
        "[Heavy ensemble score coverage in candidate pool]",
        (
            f"non_null: {int(base_cand['heavy_ensemble_score'].notna().sum())} / {len(base_cand)} "
            f"({base_cand['heavy_ensemble_score'].notna().mean():.3f})"
            if "heavy_ensemble_score" in base_cand.columns
            else "heavy_ensemble_score column missing in candidate pool"
        ),
        "",
        "[Heavy ensemble score coverage in final library]",
        (
            f"non_null: {int(lib['heavy_ensemble_score'].notna().sum())} / {len(lib)} "
            f"({lib['heavy_ensemble_score'].notna().mean():.3f})\n"
            + pd.to_numeric(lib["heavy_ensemble_score"], errors="coerce").describe().to_string()
            if "heavy_ensemble_score" in lib.columns
            else "heavy_ensemble_score column missing in final library"
        ),
        "",
        f"Saved CSV: {OUT_CSV}",
        f"Saved FASTA: {OUT_FASTA}",
        f"Saved gene/sequence-cluster summary: {OUT_DIVERSITY}",
        f"Saved uAUG source summary: {OUT_UAUG_SUMMARY_TXT}",
        f"Saved uAUG source table: {OUT_UAUG_SUMMARY_CSV}",
        f"Saved uAUG-positive rows: {OUT_UAUG_POSITIVE_CSV}",
        f"Saved uAUG=0 production validation: {OUT_UAUG0_VALIDATION}",
        f"Saved uAUG=0 dry-run summary: {OUT_UAUG0_DRY_SUMMARY}",
        f"Saved uAUG=0 dry-run shortfall: {OUT_UAUG0_SHORTFALL}",
        f"Saved uAUG=0 replacement candidates: {OUT_UAUG0_REPLACEMENTS}",
        f"Saved evidence refill audit: {OUT_REFILL_AUDIT}",
        f"Saved v1.4 selection policy QC: {OUT_SELECTION_QC}",
    ]
    OUT_QC.write_text("\n".join(q), encoding="utf-8")
    print("[SAVED]", OUT_CSV)
    print("[SAVED]", OUT_FASTA)
    print("[SAVED]", OUT_QC)
    print("[SAVED]", OUT_DIVERSITY)
    print("[SAVED]", OUT_UAUG_SUMMARY_TXT)
    print("[SAVED]", OUT_UAUG_SUMMARY_CSV)
    print("[SAVED]", OUT_UAUG_POSITIVE_CSV)
    print("[SAVED]", OUT_UAUG0_VALIDATION)
    print("[SAVED]", OUT_UAUG0_DRY_SUMMARY)
    print("[SAVED]", OUT_UAUG0_SHORTFALL)
    print("[SAVED]", OUT_UAUG0_REPLACEMENTS)
    print("[SAVED]", OUT_REFILL_AUDIT)
    print("[SAVED]", OUT_SELECTION_QC)


if __name__ == "__main__":
    main()
