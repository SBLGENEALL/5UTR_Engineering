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
for p in [OUT_CSV.parent, OUT_FASTA.parent, OUT_QC.parent, OUT_DIVERSITY.parent]:
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
    args = ap.parse_args()

    path = choose_input(args.input)
    print("[LOAD]", path)
    df = prep(pd.read_csv(path))
    cand = df[
        df["length"].between(50, 100) &
        pd.to_numeric(df["gc_content"], errors="coerce").between(0.30, 0.75) &
        (pd.to_numeric(df["uaug_count"], errors="coerce").fillna(999) <= 1) &
        df["forbidden_sites"].eq("") &
        df[SEQ].str.len().gt(0) &
        ~df[SEQ].str.contains("N", regex=False) &
        df.get("robust_public_te_rank", pd.Series(np.nan, index=df.index)).notna()
    ].copy()
    cand = cand.drop_duplicates(subset=[SEQ]).reset_index(drop=True)
    cand["cluster_diverse_evidence_score"] = evidence_score(cand)

    print(f"[CANDIDATES] {len(cand):,}")
    if len(cand) < args.n:
        print(f"[WARN] candidate count {len(cand)} < requested {args.n}; output will be smaller unless controls add rows")

    selected = []
    used_seq = set()
    cluster_counts = defaultdict(int)

    def can_take(row, max_per):
        seq = row[SEQ]
        cid = str(row.get("seq_cluster_id", seq))
        return seq not in used_seq and cluster_counts[cid] < max_per

    def take(pool, n, group, sort_cols, ascending=None, max_per_cluster=None):
        nonlocal selected, used_seq, cluster_counts
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
        if rows:
            out = pd.DataFrame(rows)
            out["library_group"] = group
            selected.append(out)
            print(f"  take {group}: {len(out)}/{n}")
        else:
            print(f"  take {group}: 0/{n}")

    # Quotas total 1950 + up to 50 reference controls.
    quotas = {
        "A_publicTE_high_confidence": 500,
        "B_TE_model_classifier_supported": 300,
        "C_protein_abundance_supported": 250,
        "D_protein_residual_supported": 250,
        "E_multiomics_consensus_high": 250,
        "F_sequence_diverse_exploratory": 200,
        "G_length_GC_uAUG_diversity": 100,
        "H_low_signal_negative_controls": 100,
    }

    base_sort = ["cluster_diverse_evidence_score", "robust_public_te_rank"]
    take(cand, quotas["A_publicTE_high_confidence"], "A_publicTE_high_confidence", ["robust_public_te_rank", "day_consensus_TE_rank" if "day_consensus_TE_rank" in cand.columns else "cluster_diverse_evidence_score"])

    model_cols = [c for c in ["heavy_ensemble_score", "automl_ensemble_score", "proteomics_enriched_score", "model_pred_rank_40_200train"] if c in cand.columns]
    if model_cols:
        cand["model_support_score"] = cand[model_cols].apply(pd.to_numeric, errors="coerce").mean(axis=1)
    else:
        cand["model_support_score"] = cand["cluster_diverse_evidence_score"]
    take(cand, quotas["B_TE_model_classifier_supported"], "B_TE_model_classifier_supported", ["model_support_score", "robust_public_te_rank"])

    if "protein_abundance_rank" in cand.columns:
        take(cand[cand["protein_abundance_rank"].notna()], quotas["C_protein_abundance_supported"], "C_protein_abundance_supported", ["protein_abundance_rank", "cluster_diverse_evidence_score"])
    if "protein_residual_rank" in cand.columns:
        take(cand[cand["protein_residual_rank"].notna()], quotas["D_protein_residual_supported"], "D_protein_residual_supported", ["protein_residual_rank", "cluster_diverse_evidence_score"])
    if "multi_omics_utr_rank" in cand.columns:
        take(cand, quotas["E_multiomics_consensus_high"], "E_multiomics_consensus_high", ["multi_omics_utr_rank", "cluster_diverse_evidence_score"])

    # exploratory: mid/high but not top-only; use bins to diversify length and GC.
    exploratory = cand[cand["robust_public_te_rank"].between(0.35, 0.90)].copy()
    exploratory["length_bin"] = pd.cut(exploratory["length"], [49, 60, 75, 90, 100], labels=["50-60", "61-75", "76-90", "91-100"])
    exploratory["gc_bin"] = pd.cut(exploratory["gc_content"], [0.299, 0.45, 0.60, 0.75], labels=["30-45", "45-60", "60-75"])
    exploratory["diversity_key"] = exploratory["length_bin"].astype(str) + ":" + exploratory["gc_bin"].astype(str) + ":uAUG" + exploratory["uaug_count"].astype(str)
    # Interleave by diversity key
    exp_rows = []
    for _, g in exploratory.sort_values("cluster_diverse_evidence_score", ascending=False).groupby("diversity_key", dropna=False):
        exp_rows.append(g.head(max(1, quotas["F_sequence_diverse_exploratory"] // max(1, exploratory["diversity_key"].nunique()))))
    exp_pool = pd.concat(exp_rows).sort_values("cluster_diverse_evidence_score", ascending=False) if exp_rows else exploratory
    take(exp_pool, quotas["F_sequence_diverse_exploratory"], "F_sequence_diverse_exploratory", ["cluster_diverse_evidence_score", "robust_public_te_rank"])

    diversity = cand.copy()
    diversity["diversity_score"] = (
        (1 - abs(pd.to_numeric(diversity["gc_content"], errors="coerce") - 0.52)) +
        (1 - abs(pd.to_numeric(diversity["length"], errors="coerce") - 75) / 50) +
        0.2 * (diversity["uaug_count"] == 0).astype(float)
    )
    take(diversity, quotas["G_length_GC_uAUG_diversity"], "G_length_GC_uAUG_diversity", ["diversity_score", "cluster_diverse_evidence_score"])

    low = cand.copy()
    low["low_score"] = 1 - pd.to_numeric(low["robust_public_te_rank"], errors="coerce").fillna(1)
    if "protein_abundance_rank" in low.columns:
        low["low_score"] += (1 - pd.to_numeric(low["protein_abundance_rank"], errors="coerce").fillna(1)) * 0.25
    if "protein_residual_rank" in low.columns:
        low["low_score"] += (1 - pd.to_numeric(low["protein_residual_rank"], errors="coerce").fillna(1)) * 0.25
    take(low[low["robust_public_te_rank"] <= 0.30], quotas["H_low_signal_negative_controls"], "H_low_signal_negative_controls", ["low_score", "gc_content"])

    lib = pd.concat(selected, ignore_index=True) if selected else cand.head(0)

    # add reference controls up to 50 if present, without enforcing cluster.
    controls = parse_reference_controls()
    if len(controls):
        controls = prep(controls)
        controls["library_group"] = "I_reference_controls"
        controls["is_reference_control"] = True
        controls["cluster_diverse_evidence_score"] = np.nan
        controls = controls.head(50)
        lib["is_reference_control"] = False
        lib = pd.concat([lib, controls], ignore_index=True, sort=False)
        print(f"  added reference controls: {len(controls)}")
    else:
        lib["is_reference_control"] = False
        print("  reference controls: none found")

    # Final fill if below requested n.
    if len(lib) < args.n:
        used_seq = set(lib[SEQ])
        cluster_counts = defaultdict(int, lib[~lib["is_reference_control"].fillna(False)].groupby("seq_cluster_id").size().to_dict() if "seq_cluster_id" in lib.columns else {})
        fill = []
        pool = cand[~cand[SEQ].isin(used_seq)].sort_values(["cluster_diverse_evidence_score", "robust_public_te_rank"], ascending=False)
        for _, row in pool.iterrows():
            if len(lib) + len(fill) >= args.n:
                break
            cid = str(row.get("seq_cluster_id", row[SEQ]))
            if cluster_counts[cid] < args.allow_cluster_fill:
                fill.append(row)
                cluster_counts[cid] += 1
        if fill:
            f = pd.DataFrame(fill)
            f["library_group"] = "J_fill_best_remaining_allow_cluster2"
            f["is_reference_control"] = False
            lib = pd.concat([lib, f], ignore_index=True, sort=False)
            print(f"  fill: {len(f)}")

    lib = lib.drop_duplicates(subset=[SEQ], keep="first").head(args.n).copy()
    lib["library_index"] = np.arange(1, len(lib) + 1)

    front = [c for c in [
        "library_index", "library_group", "is_reference_control", "utr_id", "gene_id", "gene_name", SEQ,
        "length", "gc_content", "uaug_count", "seq_cluster_id", "seq_cluster_size", "forbidden_sites",
        "robust_public_te_rank", "day_consensus_TE_rank", "protein_abundance_rank", "protein_residual_rank", "multi_omics_utr_rank",
        "cluster_diverse_evidence_score", "model_support_score", "has_proteomics_label"
    ] if c in lib.columns]
    lib = lib[front + [c for c in lib.columns if c not in front]]
    lib.to_csv(OUT_CSV, index=False)
    write_fasta(lib, OUT_FASTA)
    write_final_diversity_summary(lib, OUT_DIVERSITY)

    q = [
        "Cluster-diverse evidence-balanced 2000 library summary",
        "=" * 100,
        f"input: {path}",
        f"candidate_pool_after_QC: {len(cand)}",
        f"requested_n: {args.n}",
        f"selected_n: {len(lib)}",
        f"max_per_cluster_primary: {args.max_per_cluster}",
        f"allow_cluster_fill: {args.allow_cluster_fill}",
        "",
        "[Group counts]",
        lib["library_group"].value_counts(dropna=False).to_string(),
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
        "[Heavy ensemble score coverage in candidate pool]",
        (
            f"non_null: {int(cand['heavy_ensemble_score'].notna().sum())} / {len(cand)} "
            f"({cand['heavy_ensemble_score'].notna().mean():.3f})"
            if "heavy_ensemble_score" in cand.columns
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
    ]
    OUT_QC.write_text("\n".join(q), encoding="utf-8")
    print("[SAVED]", OUT_CSV)
    print("[SAVED]", OUT_FASTA)
    print("[SAVED]", OUT_QC)
    print("[SAVED]", OUT_DIVERSITY)


if __name__ == "__main__":
    main()
