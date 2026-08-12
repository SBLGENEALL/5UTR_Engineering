"""Sort-seq bin expression scoring.

Consumes the combined amplicon-QC output of NGS_LibraryQC (variant_count_matrix.csv,
produced by running amplicon_qc.py / amplicon_qc_parallel.py once across pDNA,
unsorted, and bin1..bin6 gDNA FASTQs against the same 5'UTR reference used for the
original plasmid-pool NGS run) plus a bin design table (FACS-recorded cells sorted
and representative MFI per bin) to compute, per library variant:

  - a cell-count-weighted expression score  = weighted mean of log10(bin MFI),
    weighted by each variant's estimated fraction of cells in each bin
    (Sort-seq / FACS-seq scoring, Sharon et al. 2012 Nat Biotechnol;
    Kosuri et al. 2013 PNAS)
  - the spread (weighted stdev) of that per-variant bin distribution, as a
    noise/bimodality indicator
  - QC flags tracing each variant through pDNA -> unsorted -> bins, so
    "enriched in a high bin" claims are only made for variants that survived
    cloning and transfection/selection with adequate read support
  - optional correlation against this project's own model-predicted rank
    (step 10 library design output), to validate whether predicted-high-TE
    variants actually land in high-MFI bins
"""
from __future__ import annotations

from pathlib import Path
import argparse
import re
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

BASE = Path.cwd()
SEQ = "utr5_sequence_tss_corrected"

DEFAULT_COUNT_MATRIX_CANDIDATES = [
    BASE / "results/combined/variant_count_matrix.csv",
    BASE / "08_sortseq_validation/input/variant_count_matrix.csv",
]
DEFAULT_BIN_METADATA = BASE / "01_pipeline/config/sortseq_bin_design.csv"
DEFAULT_LIBRARY_DESIGN = BASE / "07_library_design/tables/selected_2000_50_100bp_cluster_diverse_evidence_balanced_library.csv"

OUTDIR = BASE / "08_sortseq_validation"
TABLEDIR = OUTDIR / "tables"
QCDIR = OUTDIR / "qc"
PLOTDIR = OUTDIR / "plots"
for d in [TABLEDIR, QCDIR, PLOTDIR]:
    d.mkdir(parents=True, exist_ok=True)


def clean_seq(x):
    if pd.isna(x):
        return ""
    return re.sub(r"[^ACGTN]", "", str(x).upper().replace("U", "T"))


def find_file(path_arg, candidates, label):
    if path_arg:
        p = Path(path_arg)
        if not p.exists():
            raise SystemExit(f"{label} not found: {p}")
        return p
    for p in candidates:
        if p.exists():
            return p
    raise SystemExit(f"No {label} found. Pass an explicit path.")


def load_count_matrix(path: Path) -> tuple[pd.DataFrame, list[str]]:
    df = pd.read_csv(path)
    meta_cols = [c for c in ["sequence_key", "variant_ids", "target_sequence", "length", "gc_percent"] if c in df.columns]
    sample_cols = [c for c in df.columns if c not in meta_cols]
    if not sample_cols:
        raise SystemExit(f"No sample columns found in {path}")
    for c in sample_cols:
        df[c] = pd.to_numeric(df[c], errors="coerce").fillna(0.0)
    return df, sample_cols


def load_bin_metadata(path: Path, sample_cols: list[str]) -> pd.DataFrame:
    meta = pd.read_csv(path)
    required = {"sample", "role"}
    missing = required - set(meta.columns)
    if missing:
        raise SystemExit(f"bin metadata missing columns: {missing}")
    meta["role"] = meta["role"].str.strip().str.lower()
    unknown = set(meta["sample"]) - set(sample_cols)
    if unknown:
        print(f"[WARN] bin metadata samples not present in count matrix columns, dropping: {sorted(unknown)}")
        meta = meta[meta["sample"].isin(sample_cols)].copy()
    bins = meta[meta["role"] == "bin"].copy()
    if bins.empty:
        raise SystemExit("bin metadata has no rows with role == 'bin'")
    for col in ["bin_order", "n_cells_sorted", "bin_mfi_geomean"]:
        bins[col] = pd.to_numeric(bins[col], errors="coerce")
    if bins["n_cells_sorted"].isna().any():
        bad = bins.loc[bins["n_cells_sorted"].isna(), "sample"].tolist()
        raise SystemExit(f"n_cells_sorted missing for bins: {bad}. Fill this in from the FACS sort report before scoring.")
    if bins["bin_mfi_geomean"].isna().any():
        bad = bins.loc[bins["bin_mfi_geomean"].isna(), "sample"].tolist()
        raise SystemExit(f"bin_mfi_geomean missing for bins: {bad}. Fill this in from the FACS analysis (geometric mean fluorescence per gate) before scoring.")
    bins = bins.sort_values("bin_order").reset_index(drop=True)
    pdna_samples = meta.loc[meta["role"] == "pdna", "sample"].tolist()
    unsorted_samples = meta.loc[meta["role"] == "unsorted", "sample"].tolist()
    return bins, pdna_samples, unsorted_samples


def score_variants(
    df: pd.DataFrame,
    sample_totals: dict[str, float],
    bins: pd.DataFrame,
    pdna_samples: list[str],
    unsorted_samples: list[str],
    min_reads_pdna: int,
    min_reads_unsorted: int,
    min_bin_reads: int,
) -> pd.DataFrame:
    bin_samples = bins["sample"].tolist()
    bin_mfi = dict(zip(bins["sample"], bins["bin_mfi_geomean"]))
    bin_cells = dict(zip(bins["sample"], bins["n_cells_sorted"]))
    log_mfi = {s: np.log10(bin_mfi[s]) for s in bin_samples}

    out = df.copy()

    out["pdna_reads"] = out[pdna_samples].sum(axis=1) if pdna_samples else np.nan
    out["unsorted_reads"] = out[unsorted_samples].sum(axis=1) if unsorted_samples else np.nan
    out["pdna_rpm"] = 1e6 * out["pdna_reads"] / sum(sample_totals[s] for s in pdna_samples) if pdna_samples else np.nan
    out["unsorted_rpm"] = 1e6 * out["unsorted_reads"] / sum(sample_totals[s] for s in unsorted_samples) if unsorted_samples else np.nan

    # estimated cells of each variant per bin = (variant read fraction in that
    # bin's library) x (total cells actually sorted into that bin)
    est_cells = pd.DataFrame(index=out.index)
    for s in bin_samples:
        read_frac = out[s] / sample_totals[s] if sample_totals[s] else 0.0
        est_cells[s] = read_frac * bin_cells[s]

    out["bin_reads_total"] = out[bin_samples].sum(axis=1)
    row_total = est_cells.sum(axis=1)
    has_signal = row_total > 0

    p = est_cells.div(row_total.replace(0, np.nan), axis=0)
    log_mfi_vec = np.array([log_mfi[s] for s in bin_samples])
    score = (p[bin_samples].values * log_mfi_vec).sum(axis=1)
    mean_score = np.where(has_signal, score, np.nan)
    diff_sq = (p[bin_samples].values - 0) * 0  # placeholder to keep shape
    spread = np.sqrt(
        np.nansum(p[bin_samples].values * (log_mfi_vec - mean_score.reshape(-1, 1)) ** 2, axis=1)
    )
    dominant_bin = p[bin_samples].idxmax(axis=1)

    out["sortseq_log10_mfi_score"] = mean_score
    out["sortseq_score_spread"] = np.where(has_signal, spread, np.nan)
    out["sortseq_dominant_bin"] = np.where(has_signal, dominant_bin, None)
    for s in bin_samples:
        out[f"bin_fraction__{s}"] = p[s]

    if pdna_samples:
        out["flag_dropped_at_cloning"] = out["pdna_reads"] < min_reads_pdna
    else:
        out["flag_dropped_at_cloning"] = False

    if pdna_samples and unsorted_samples:
        out["flag_lost_at_transfection_or_selection"] = (
            (out["pdna_reads"] >= min_reads_pdna) & (out["unsorted_reads"] < min_reads_unsorted)
        )
    else:
        out["flag_lost_at_transfection_or_selection"] = False

    out["flag_low_bin_coverage"] = out["bin_reads_total"] < min_bin_reads
    out["flag_retained_high_confidence"] = (
        (~out["flag_dropped_at_cloning"].astype(bool))
        & (~out["flag_lost_at_transfection_or_selection"].astype(bool))
        & (~out["flag_low_bin_coverage"].astype(bool))
        & has_signal
    )
    return out


def attach_design_rank(scored: pd.DataFrame, design_path: Path | None) -> pd.DataFrame:
    if design_path is None or not Path(design_path).exists():
        print("[INFO] no library-design table given/found; skipping predicted-rank correlation")
        return scored
    design = pd.read_csv(design_path)
    if SEQ not in design.columns or "target_sequence" not in scored.columns:
        print(f"[WARN] cannot join design table: need '{SEQ}' in design and 'target_sequence' in count matrix")
        return scored
    design = design.copy()
    design["_seq_key"] = design[SEQ].map(clean_seq)
    scored = scored.copy()
    scored["_seq_key"] = scored["target_sequence"].map(clean_seq)
    keep_cols = [c for c in [
        "_seq_key", "gene_name", "library_group", "library_index",
        "robust_public_te_rank", "cluster_diverse_evidence_score", "heavy_ensemble_score",
    ] if c in design.columns]
    merged = scored.merge(design[keep_cols], on="_seq_key", how="left")
    return merged.drop(columns=["_seq_key"])


def spearman(a: pd.Series, b: pd.Series) -> tuple[float, int]:
    d = pd.DataFrame({"a": a, "b": b}).dropna()
    if len(d) < 3:
        return float("nan"), len(d)
    return d["a"].corr(d["b"], method="spearman"), len(d)


def write_plots(scored: pd.DataFrame, bins: pd.DataFrame):
    bin_samples = bins["sample"].tolist()

    # 1) predicted rank vs measured sort-seq score
    if "robust_public_te_rank" in scored.columns:
        fig, ax = plt.subplots(figsize=(6, 6))
        d = scored[scored["flag_retained_high_confidence"]]
        ax.scatter(d["robust_public_te_rank"], d["sortseq_log10_mfi_score"], s=8, alpha=0.4)
        r, n = spearman(d["robust_public_te_rank"], d["sortseq_log10_mfi_score"])
        ax.set_xlabel("model-predicted TE rank (step 10 input)")
        ax.set_ylabel("measured Sort-seq score (weighted log10 MFI)")
        ax.set_title(f"Predicted rank vs measured expression (Spearman r={r:.3f}, n={n})")
        fig.tight_layout()
        fig.savefig(PLOTDIR / "01_predicted_rank_vs_measured_score.png", dpi=150)
        plt.close(fig)

    # 2) pDNA vs unsorted representation retention
    if "pdna_rpm" in scored.columns and "unsorted_rpm" in scored.columns:
        fig, ax = plt.subplots(figsize=(6, 6))
        d = scored[(scored["pdna_rpm"] > 0) & (scored["unsorted_rpm"] > 0)]
        ax.scatter(np.log10(d["pdna_rpm"]), np.log10(d["unsorted_rpm"]), s=8, alpha=0.4)
        r, n = spearman(d["pdna_rpm"], d["unsorted_rpm"])
        ax.set_xlabel("log10(pDNA RPM)")
        ax.set_ylabel("log10(unsorted gDNA RPM)")
        ax.set_title(f"Library representation retained through transfection (Spearman r={r:.3f}, n={n})")
        fig.tight_layout()
        fig.savefig(PLOTDIR / "02_pdna_vs_unsorted_representation.png", dpi=150)
        plt.close(fig)

    # 3) stacked bin-fraction bars for top/bottom 20 variants by score
    d = scored.dropna(subset=["sortseq_log10_mfi_score"]).sort_values("sortseq_log10_mfi_score", ascending=False)
    top = d.head(20)
    bottom = d.tail(20)
    for name, sub in [("top20_highest_score", top), ("bottom20_lowest_score", bottom)]:
        fig, ax = plt.subplots(figsize=(10, 5))
        bottom_stack = np.zeros(len(sub))
        labels = sub.get("variant_ids", sub.get("sequence_key", pd.Series(range(len(sub))))).astype(str).str.slice(0, 20)
        for s in bin_samples:
            vals = sub[f"bin_fraction__{s}"].fillna(0).values
            ax.bar(range(len(sub)), vals, bottom=bottom_stack, label=s)
            bottom_stack += vals
        ax.set_xticks(range(len(sub)))
        ax.set_xticklabels(labels, rotation=90, fontsize=6)
        ax.set_ylabel("estimated fraction of variant's cells per bin")
        ax.set_title(f"Bin distribution: {name}")
        ax.legend(fontsize=6, ncol=3)
        fig.tight_layout()
        fig.savefig(PLOTDIR / f"03_bin_distribution_{name}.png", dpi=150)
        plt.close(fig)

    # 4) score distribution histogram
    fig, ax = plt.subplots(figsize=(6, 4))
    ax.hist(scored["sortseq_log10_mfi_score"].dropna(), bins=40)
    ax.set_xlabel("Sort-seq score (weighted log10 MFI)")
    ax.set_ylabel("variant count")
    ax.set_title("Distribution of measured expression scores")
    fig.tight_layout()
    fig.savefig(PLOTDIR / "04_score_distribution.png", dpi=150)
    plt.close(fig)


def main():
    ap = argparse.ArgumentParser(description="Score 5'UTR library variants from Sort-seq bin NGS (pDNA/unsorted/bin1..N)")
    ap.add_argument("--count-matrix", default=None, help="NGS_LibraryQC combined/variant_count_matrix.csv (raw counts, not RPM)")
    ap.add_argument("--bin-metadata", default=None, help="CSV: sample,role,bin_order,gate_pct,n_cells_sorted,bin_mfi_geomean")
    ap.add_argument("--library-design", default=None, help="step-10 selected library CSV, for predicted-rank correlation")
    ap.add_argument("--min-reads-pdna", type=int, default=20)
    ap.add_argument("--min-reads-unsorted", type=int, default=20)
    ap.add_argument("--min-bin-reads", type=int, default=20, help="minimum total reads across all bins to trust a variant's score")
    ap.add_argument("--library-size", type=int, default=None, help="intended library complexity, for a retrospective coverage check")
    args = ap.parse_args()

    count_path = find_file(args.count_matrix, DEFAULT_COUNT_MATRIX_CANDIDATES, "variant_count_matrix.csv")
    bin_meta_path = find_file(args.bin_metadata, [DEFAULT_BIN_METADATA], "bin metadata CSV")
    print("[LOAD]", count_path)
    df, sample_cols = load_count_matrix(count_path)
    sample_totals = {c: df[c].sum() for c in sample_cols}

    print("[LOAD]", bin_meta_path)
    bins, pdna_samples, unsorted_samples = load_bin_metadata(bin_meta_path, sample_cols)
    print(f"[DESIGN] pDNA samples: {pdna_samples} | unsorted samples: {unsorted_samples} | bins (low->high MFI): {bins['sample'].tolist()}")

    if args.library_size:
        for _, row in bins.iterrows():
            per_variant = row["n_cells_sorted"] / args.library_size
            flag = "OK" if per_variant >= 100 else "LOW COVERAGE"
            print(f"[COVERAGE] {row['sample']}: {row['n_cells_sorted']:,.0f} cells / {args.library_size} variants = {per_variant:,.0f}x  [{flag}]")

    scored = score_variants(
        df, sample_totals, bins, pdna_samples, unsorted_samples,
        args.min_reads_pdna, args.min_reads_unsorted, args.min_bin_reads,
    )

    design_path = Path(args.library_design) if args.library_design else DEFAULT_LIBRARY_DESIGN
    scored = attach_design_rank(scored, design_path)

    out_csv = TABLEDIR / "variant_sortseq_expression_score.csv"
    scored.to_csv(out_csv, index=False)
    print("[WRITE]", out_csv, f"({len(scored):,} variants)")

    n_total = len(scored)
    n_dropped = int(scored["flag_dropped_at_cloning"].sum()) if pdna_samples else None
    n_lost = int(scored["flag_lost_at_transfection_or_selection"].sum())
    n_low_cov = int(scored["flag_low_bin_coverage"].sum())
    n_retained = int(scored["flag_retained_high_confidence"].sum())
    r_pdna_unsorted, n_pu = (spearman(scored["pdna_rpm"], scored["unsorted_rpm"]) if pdna_samples and unsorted_samples else (float("nan"), 0))
    r_rank_score, n_rs = (spearman(scored["robust_public_te_rank"], scored["sortseq_log10_mfi_score"])
                           if "robust_public_te_rank" in scored.columns else (float("nan"), 0))

    summary_lines = [
        f"total_variants={n_total}",
        f"dropped_at_cloning(pDNA<{args.min_reads_pdna}reads)={n_dropped}",
        f"lost_at_transfection_or_selection(unsorted<{args.min_reads_unsorted}reads)={n_lost}",
        f"low_bin_coverage(<{args.min_bin_reads}total bin reads)={n_low_cov}",
        f"retained_high_confidence={n_retained} ({100*n_retained/n_total:.1f}%)",
        f"spearman_pDNA_vs_unsorted_rpm={r_pdna_unsorted:.3f} (n={n_pu})",
        f"spearman_predicted_rank_vs_measured_score={r_rank_score:.3f} (n={n_rs})",
    ]
    qc_path = QCDIR / "sortseq_qc_summary.txt"
    qc_path.write_text("\n".join(summary_lines) + "\n", encoding="utf-8")
    print("[WRITE]", qc_path)
    for line in summary_lines:
        print(" ", line)

    write_plots(scored, bins)
    print("[DONE] plots ->", PLOTDIR)


if __name__ == "__main__":
    main()
