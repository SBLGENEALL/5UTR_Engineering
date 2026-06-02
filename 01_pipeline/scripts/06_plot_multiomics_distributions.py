from pathlib import Path
import argparse
import re
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

BASE = Path.cwd()
SEQ = "utr5_sequence_tss_corrected"

DEFAULT_INPUTS = [
    BASE / "04_te_labeling/tables/tss_corrected_5utr_multiomics_labels.csv",
    BASE / "04_te_labeling/tables/tss_corrected_5utr_robust_public_te_labels.csv",
]

PLOTDIR = BASE / "06_modeling/plots/multiomics_distributions"
QCDIR = BASE / "04_te_labeling/qc"
TABLEDIR = BASE / "06_modeling/tables"
for d in [PLOTDIR, QCDIR, TABLEDIR]:
    d.mkdir(parents=True, exist_ok=True)


def clean_seq(x):
    if pd.isna(x):
        return ""
    return re.sub(r"[^ACGTN]", "", str(x).upper().replace("U", "T"))


def gc_content(seq):
    seq = clean_seq(seq)
    return (seq.count("G") + seq.count("C")) / len(seq) if len(seq) else np.nan


def find_input(path_arg=None):
    if path_arg:
        p = Path(path_arg)
        if not p.exists():
            raise SystemExit(f"Input file not found: {p}")
        return p
    for p in DEFAULT_INPUTS:
        if p.exists():
            return p
    raise SystemExit("No multiomics/publicTE label file found.")


def numeric(df, col):
    return pd.to_numeric(df[col], errors="coerce")


def add_basic_columns(df):
    df = df.copy()
    if SEQ in df.columns:
        df[SEQ] = df[SEQ].map(clean_seq)

    if "utr5_length_final" in df.columns:
        df["length_for_plot"] = pd.to_numeric(df["utr5_length_final"], errors="coerce")
    elif "utr5_length_tss_corrected" in df.columns:
        df["length_for_plot"] = pd.to_numeric(df["utr5_length_tss_corrected"], errors="coerce")
    elif SEQ in df.columns:
        df["length_for_plot"] = df[SEQ].str.len()
    else:
        df["length_for_plot"] = np.nan

    if "gc_content" not in df.columns and SEQ in df.columns:
        df["gc_content"] = df[SEQ].map(gc_content)

    if "uaug_count" not in df.columns and SEQ in df.columns:
        df["uaug_count"] = df[SEQ].str.count("ATG")

    if "has_proteomics_label" not in df.columns:
        prot_cols = [c for c in ["protein_abundance_rank", "protein_residual_rank", "abundance_proxy"] if c in df.columns]
        df["has_proteomics_label"] = df[prot_cols].notna().any(axis=1) if prot_cols else False

    return df


def existing(df, cols):
    return [c for c in cols if c in df.columns]


def plot_hist(df, col, suffix=""):
    x = numeric(df, col).dropna()
    if len(x) == 0:
        return
    plt.figure(figsize=(7, 5))
    plt.hist(x, bins=50)
    plt.xlabel(col)
    plt.ylabel("Count")
    plt.title(f"Distribution: {col}{suffix}\nn={len(x):,}")
    plt.tight_layout()
    plt.savefig(PLOTDIR / f"hist_{col}{suffix.replace(' ', '_')}.png", dpi=220)
    plt.close()


def plot_ecdf(df, col, suffix=""):
    x = np.sort(numeric(df, col).dropna().values)
    if len(x) == 0:
        return
    y = np.arange(1, len(x) + 1) / len(x)
    plt.figure(figsize=(7, 5))
    plt.plot(x, y)
    plt.xlabel(col)
    plt.ylabel("Cumulative fraction")
    plt.title(f"ECDF: {col}{suffix}\nn={len(x):,}")
    plt.tight_layout()
    plt.savefig(PLOTDIR / f"ecdf_{col}{suffix.replace(' ', '_')}.png", dpi=220)
    plt.close()


def plot_overlay_hist(df, cols, filename, title):
    cols = [c for c in cols if c in df.columns and numeric(df, c).notna().sum() > 0]
    if not cols:
        return
    plt.figure(figsize=(8, 5))
    for c in cols:
        x = numeric(df, c).dropna()
        plt.hist(x, bins=40, alpha=0.35, label=f"{c} n={len(x):,}")
    plt.xlabel("Rank / score value")
    plt.ylabel("Count")
    plt.title(title)
    plt.legend(fontsize=8)
    plt.tight_layout()
    plt.savefig(PLOTDIR / filename, dpi=220)
    plt.close()


def plot_box(df, cols, filename, title):
    data, labels = [], []
    for c in cols:
        if c in df.columns:
            x = numeric(df, c).dropna().values
            if len(x):
                data.append(x)
                labels.append(c)
    if not data:
        return
    plt.figure(figsize=(max(8, 0.7 * len(labels)), 5))
    plt.boxplot(data, labels=labels, showfliers=False)
    plt.ylabel("Rank / score")
    plt.title(title)
    plt.xticks(rotation=45, ha="right")
    plt.tight_layout()
    plt.savefig(PLOTDIR / filename, dpi=220)
    plt.close()


def plot_corr_heatmap(df, cols, filename, title):
    cols = [c for c in cols if c in df.columns and numeric(df, c).notna().sum() >= 20]
    if len(cols) < 2:
        return
    tmp = df[cols].apply(pd.to_numeric, errors="coerce")
    corr = tmp.corr(method="spearman")
    corr.to_csv(TABLEDIR / filename.replace(".png", ".csv"))
    plt.figure(figsize=(max(7, 0.7 * len(cols)), max(6, 0.7 * len(cols))))
    im = plt.imshow(corr.values, vmin=-1, vmax=1)
    plt.colorbar(im, fraction=0.046, pad=0.04, label="Spearman correlation")
    plt.xticks(range(len(cols)), cols, rotation=45, ha="right")
    plt.yticks(range(len(cols)), cols)
    plt.title(title)
    for i in range(len(cols)):
        for j in range(len(cols)):
            plt.text(j, i, f"{corr.iloc[i, j]:.2f}", ha="center", va="center", fontsize=8)
    plt.tight_layout()
    plt.savefig(PLOTDIR / filename, dpi=220)
    plt.close()


def plot_scatter(df, xcol, ycol, sample_n=15000):
    if xcol not in df.columns or ycol not in df.columns:
        return
    tmp = df[[xcol, ycol]].apply(pd.to_numeric, errors="coerce").dropna()
    if len(tmp) < 20:
        return
    if len(tmp) > sample_n:
        tmp = tmp.sample(sample_n, random_state=42)
    corr = tmp.corr(method="spearman").iloc[0, 1]
    plt.figure(figsize=(6, 6))
    plt.scatter(tmp[xcol], tmp[ycol], s=8, alpha=0.35)
    plt.xlabel(xcol)
    plt.ylabel(ycol)
    plt.title(f"{xcol} vs {ycol}\nSpearman={corr:.3f}, n={len(tmp):,}")
    plt.tight_layout()
    plt.savefig(PLOTDIR / f"scatter_{xcol}_vs_{ycol}.png", dpi=220)
    plt.close()


def plot_missingness(df, cols):
    rows = []
    for c in cols:
        if c in df.columns:
            rows.append((c, int(df[c].notna().sum()), int(df[c].isna().sum())))
    if not rows:
        return
    s = pd.DataFrame(rows, columns=["column", "non_missing", "missing"])
    s.to_csv(TABLEDIR / "multiomics_label_missingness.csv", index=False)
    plt.figure(figsize=(max(8, 0.45 * len(s)), 5))
    plt.bar(s["column"], s["non_missing"])
    plt.ylabel("Non-missing rows")
    plt.title("Non-missing count by label column")
    plt.xticks(rotation=45, ha="right")
    plt.tight_layout()
    plt.savefig(PLOTDIR / "missingness_non_missing_counts.png", dpi=220)
    plt.close()


def plot_proteomics_coverage_by_length(df):
    if "length_for_plot" not in df.columns or "has_proteomics_label" not in df.columns:
        return
    x = df[["length_for_plot", "has_proteomics_label"]].copy()
    x["length_for_plot"] = pd.to_numeric(x["length_for_plot"], errors="coerce")
    x = x.dropna(subset=["length_for_plot"])
    if len(x) == 0:
        return
    bins = [0, 20, 40, 50, 75, 100, 150, 200, 300, 500, 1000, np.inf]
    labels = ["0-20", "20-40", "40-50", "50-75", "75-100", "100-150", "150-200", "200-300", "300-500", "500-1000", ">1000"]
    x["length_bin"] = pd.cut(x["length_for_plot"], bins=bins, labels=labels, right=False)
    tab = x.groupby("length_bin", observed=False)["has_proteomics_label"].agg(["count", "sum"])
    tab["proteomics_fraction"] = tab["sum"] / tab["count"]
    tab.to_csv(TABLEDIR / "multiomics_proteomics_coverage_by_length.csv")
    plt.figure(figsize=(9, 5))
    plt.bar(tab.index.astype(str), tab["proteomics_fraction"])
    plt.ylabel("Fraction with proteomics label")
    plt.xlabel("UTR length bin")
    plt.title("Proteomics label coverage by UTR length")
    plt.xticks(rotation=45, ha="right")
    plt.tight_layout()
    plt.savefig(PLOTDIR / "proteomics_coverage_by_length.png", dpi=220)
    plt.close()


def quantile_summary(df, cols, name):
    rows = []
    for c in cols:
        if c not in df.columns:
            continue
        x = numeric(df, c)
        valid = x.notna()
        if valid.sum() == 0:
            continue
        rows.append({
            "column": c,
            "n_non_missing": int(valid.sum()),
            "missing": int((~valid).sum()),
            "mean": float(x.mean()),
            "std": float(x.std()),
            "min": float(x.min()),
            "q05": float(x.quantile(0.05)),
            "q20": float(x.quantile(0.20)),
            "q33": float(x.quantile(0.333)),
            "median": float(x.quantile(0.50)),
            "q67": float(x.quantile(0.667)),
            "q80": float(x.quantile(0.80)),
            "q95": float(x.quantile(0.95)),
            "max": float(x.max()),
            "bottom20_n": int((x <= x.quantile(0.20)).sum()),
            "middle60_n": int(((x > x.quantile(0.20)) & (x < x.quantile(0.80))).sum()),
            "top20_n": int((x >= x.quantile(0.80)).sum()),
            "tertile_low_n": int((x <= x.quantile(0.333)).sum()),
            "tertile_mid_n": int(((x > x.quantile(0.333)) & (x < x.quantile(0.667))).sum()),
            "tertile_high_n": int((x >= x.quantile(0.667)).sum()),
        })
    out = pd.DataFrame(rows)
    out.to_csv(TABLEDIR / name, index=False)
    return out


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--input", default=None)
    ap.add_argument("--sample-n", type=int, default=15000)
    args = ap.parse_args()

    path = find_input(args.input)
    print("[LOAD]", path)
    df = add_basic_columns(pd.read_csv(path))

    rank_cols = existing(df, [
        "robust_public_te_rank", "mean_TE_rank", "day3_TE_rank", "day6_TE_rank",
        "day_consensus_TE_rank", "residual_TE_rank", "ribo_abundance_rank",
        "protein_abundance_rank", "protein_residual_rank", "multi_omics_utr_rank",
        "proteomics_enriched_score", "automl_ensemble_score", "heavy_ensemble_score",
    ])
    raw_cols = existing(df, [
        "log2_te_mean_norm", "log2_te_day3_norm", "log2_te_day6_norm",
        "te_mean_norm", "te_day3_norm", "te_day6_norm",
        "abundance_proxy", "protein_abundance_proxy",
        "rna_day3_cpm_mean", "rna_day6_cpm_mean",
        "ribo_day3_cpm_mean", "ribo_day6_cpm_mean",
    ])
    all_cols = rank_cols + raw_cols
    key_cols = existing(df, [
        "robust_public_te_rank", "multi_omics_utr_rank", "protein_abundance_rank",
        "protein_residual_rank", "residual_TE_rank", "day_consensus_TE_rank",
    ])

    print("[rank cols]", rank_cols)
    print("[raw cols]", raw_cols)

    for c in all_cols:
        plot_hist(df, c)
        plot_ecdf(df, c)

    plot_overlay_hist(df, key_cols, "overlay_hist_key_rank_columns.png", "Key multi-omics rank distributions")
    plot_box(df, key_cols, "boxplot_key_rank_columns.png", "Key multi-omics rank distributions")
    plot_corr_heatmap(df, key_cols, "spearman_correlation_heatmap_key_rank_columns.png", "Spearman correlation among key ranks")
    plot_missingness(df, all_cols)
    plot_proteomics_coverage_by_length(df)

    for xcol, ycol in [
        ("robust_public_te_rank", "multi_omics_utr_rank"),
        ("robust_public_te_rank", "protein_abundance_rank"),
        ("robust_public_te_rank", "protein_residual_rank"),
        ("protein_abundance_rank", "protein_residual_rank"),
        ("residual_TE_rank", "protein_residual_rank"),
        ("day_consensus_TE_rank", "protein_abundance_rank"),
        ("log2_te_mean_norm", "protein_abundance_rank"),
        ("log2_te_mean_norm", "protein_residual_rank"),
    ]:
        plot_scatter(df, xcol, ycol, sample_n=args.sample_n)

    mask = pd.Series(True, index=df.index)
    if "length_for_plot" in df.columns:
        mask &= pd.to_numeric(df["length_for_plot"], errors="coerce").between(50, 100)
    if "gc_content" in df.columns:
        mask &= pd.to_numeric(df["gc_content"], errors="coerce").between(0.30, 0.75)
    if "uaug_count" in df.columns:
        mask &= pd.to_numeric(df["uaug_count"], errors="coerce").fillna(999) <= 1
    subset = df[mask].copy()

    if len(subset):
        plot_overlay_hist(subset, key_cols, "overlay_hist_key_rank_columns_50_100_QC_subset.png", f"Key ranks in 50-100 bp QC subset, n={len(subset):,}")
        plot_box(subset, key_cols, "boxplot_key_rank_columns_50_100_QC_subset.png", f"Key ranks in 50-100 bp QC subset, n={len(subset):,}")
        plot_corr_heatmap(subset, key_cols, "spearman_correlation_heatmap_50_100_QC_subset.png", f"Spearman correlation in 50-100 bp QC subset, n={len(subset):,}")

    q_all = quantile_summary(df, all_cols, "multiomics_distribution_quantile_summary.csv")
    q_sub = quantile_summary(subset, all_cols, "multiomics_distribution_quantile_summary_50_100_QC_subset.csv") if len(subset) else pd.DataFrame()

    summary = [
        "Multi-omics distribution QC summary",
        "=" * 90,
        f"input_file: {path}",
        f"total_rows: {len(df)}",
        f"50_100_QC_subset_rows: {len(subset)}",
        "",
        "[Columns plotted]",
        "rank_cols: " + ", ".join(rank_cols),
        "raw_cols: " + ", ".join(raw_cols),
        "",
        "[Proteomics label coverage]",
        str(df["has_proteomics_label"].value_counts(dropna=False)) if "has_proteomics_label" in df.columns else "NA",
        "",
        "[Quantile summary: all rows]",
        q_all.to_string(index=False) if len(q_all) else "No numeric label columns.",
        "",
        "[Important interpretation]",
        "- Rank columns often look uniform because they are percentile ranks.",
        "- Raw TE values may be skewed or concentrated.",
        "- If robust_public_te_rank and multi_omics_utr_rank correlation is very high, multi-omics is mostly TE-driven.",
        "- If proteomics label coverage is low in the 50-100 bp QC subset, proteomics should be used as quota/orthogonal prior.",
        "",
        f"Plots saved to: {PLOTDIR}",
        f"Tables saved to: {TABLEDIR}",
    ]
    summary_path = QCDIR / "multiomics_distribution_qc_summary.txt"
    summary_path.write_text("\n".join(summary), encoding="utf-8")

    print("[SAVED]", summary_path)
    print("[SAVED plots]", PLOTDIR)
    print("[SAVED tables]", TABLEDIR)


if __name__ == "__main__":
    main()
