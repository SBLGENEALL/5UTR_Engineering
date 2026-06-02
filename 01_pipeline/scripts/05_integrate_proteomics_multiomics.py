from __future__ import annotations

from pathlib import Path
import argparse
import json
import re
import numpy as np
import pandas as pd

BASE = Path.cwd()

LABEL_PATH = BASE / "04_te_labeling/tables/tss_corrected_5utr_robust_public_te_labels.csv"
CONFIG_PATH = BASE / "01_pipeline/config/proteomics_config.json"

OUT_LABEL = BASE / "04_te_labeling/tables/tss_corrected_5utr_multiomics_labels.csv"
OUT_READY_40_200 = BASE / "04_te_labeling/tables/tss_corrected_5utr_multiomics_training_ready_40_200bp.csv"
OUT_READY_50_100 = BASE / "04_te_labeling/tables/tss_corrected_5utr_multiomics_selection_ready_50_100bp.csv"
OUT_PROT_SUMMARY = BASE / "04_te_labeling/qc/proteomics_mapping_summary.txt"

OUT_LABEL.parent.mkdir(parents=True, exist_ok=True)
OUT_PROT_SUMMARY.parent.mkdir(parents=True, exist_ok=True)

def norm_symbol(x):
    if pd.isna(x):
        return ""
    s = str(x).strip()
    s = re.split(r"[;,|]", s)[0].strip()
    return s.upper()

def gene_id_key(x):
    if pd.isna(x):
        return ""
    s = str(x).strip()
    # NCBI GFF often has gene-100123456, GeneID:100123456, or similar.
    nums = re.findall(r"\b\d{4,}\b", s)
    if nums:
        return nums[0]
    return s.replace("gene-", "").replace("GeneID:", "").upper()

def read_table(path, sheet_name=0):
    path = Path(path)
    if path.suffix.lower() in [".xlsx", ".xls"]:
        return pd.read_excel(path, sheet_name=sheet_name)
    for sep in [",", "\t", None]:
        for enc in ["utf-8", "utf-8-sig", "latin1", "cp949"]:
            try:
                df = pd.read_csv(path, sep=sep, engine="python", encoding=enc)
                if df.shape[1] == 1 and "," in str(df.columns[0]):
                    continue
                return df
            except Exception:
                pass
    raise RuntimeError(f"Could not read table: {path}")

def rank_pct(x):
    return pd.to_numeric(x, errors="coerce").rank(pct=True)

def residual_rank(log_rna, log_protein):
    valid = np.isfinite(log_rna) & np.isfinite(log_protein)
    res = np.full(len(log_rna), np.nan)
    if valid.sum() > 10:
        slope, intercept = np.polyfit(log_rna[valid], log_protein[valid], 1)
        res[valid] = log_protein[valid] - (slope * log_rna[valid] + intercept)
    return pd.Series(res).rank(pct=True).values

def transform_abundance(x, mode):
    x = pd.to_numeric(x, errors="coerce")
    if mode == "log1p":
        return np.log1p(x.clip(lower=0))
    if mode == "log2p1":
        return np.log2(x.clip(lower=0) + 1)
    if mode == "none":
        return x
    raise ValueError(mode)

def inspect():
    cfg = json.loads(CONFIG_PATH.read_text(encoding="utf-8"))
    p = BASE / cfg["proteomics_path"]
    df = read_table(p, cfg.get("sheet_name", 0))
    print("Proteomics file:", p)
    print("Shape:", df.shape)
    print("\nColumns:")
    for c in df.columns:
        print(" ", c)
    print("\nHead:")
    print(df.head().to_string())
    print("\nConfig gene_symbol_column:", cfg.get("gene_symbol_column"))
    print("Config gene_id_column:", cfg.get("gene_id_column"))
    print("Config abundance_columns:", cfg.get("abundance_columns"))

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--inspect", action="store_true")
    args = ap.parse_args()
    if args.inspect:
        inspect()
        return

    if not CONFIG_PATH.exists():
        raise SystemExit("Missing 01_pipeline/config/proteomics_config.json")
    cfg = json.loads(CONFIG_PATH.read_text(encoding="utf-8"))
    prot_path = BASE / cfg["proteomics_path"]

    if not LABEL_PATH.exists():
        raise SystemExit(f"Missing label file: {LABEL_PATH}")
    if not prot_path.exists():
        raise SystemExit(f"Missing proteomics file: {prot_path}")

    print("[1] Load UTR/Ribo/RNA robust public TE labels")
    lab = pd.read_csv(LABEL_PATH)
    lab["gene_symbol_key"] = lab["gene_name"].map(norm_symbol)
    lab["gene_id_key"] = lab["gene_id"].map(gene_id_key) if "gene_id" in lab.columns else ""

    print("[2] Load proteomics table:", prot_path)
    prot = read_table(prot_path, cfg.get("sheet_name", 0))
    print("  proteomics shape:", prot.shape)

    gene_col = cfg.get("gene_symbol_column", "gene_symbol_key")
    gene_id_col = cfg.get("gene_id_column", None)
    ab_cols = cfg.get("abundance_columns", ["abundance_proxy"])
    if isinstance(ab_cols, str):
        ab_cols = [ab_cols]

    if gene_col not in prot.columns and (gene_id_col is None or gene_id_col not in prot.columns):
        raise KeyError(f"Need gene_symbol_column or gene_id_column in proteomics table. Available={list(prot.columns)}")
    ab_cols = [c for c in ab_cols if c in prot.columns]
    if not ab_cols:
        raise KeyError(f"No abundance columns found from config. Available={list(prot.columns)}")

    if gene_col in prot.columns:
        prot["gene_symbol_key"] = prot[gene_col].map(norm_symbol)
    else:
        prot["gene_symbol_key"] = ""
    if gene_id_col and gene_id_col in prot.columns:
        prot["gene_id_key"] = prot[gene_id_col].map(gene_id_key)
    else:
        prot["gene_id_key"] = ""

    transformed = []
    for c in ab_cols:
        transformed.append(transform_abundance(prot[c], cfg.get("protein_abundance_transform", "none")).rename(c))
    ab = pd.concat(transformed, axis=1)
    if cfg.get("protein_aggregation", "mean") == "median":
        prot["protein_abundance_value"] = ab.median(axis=1, skipna=True)
    else:
        prot["protein_abundance_value"] = ab.mean(axis=1, skipna=True)

    prot = prot.dropna(subset=["protein_abundance_value"]).copy()

    # aggregate by gene_id if available, otherwise by symbol
    if prot["gene_id_key"].astype(str).str.len().gt(0).any():
        prot_summary_id = prot[prot["gene_id_key"].astype(str).str.len() > 0][["gene_id_key","protein_abundance_value"]].groupby("gene_id_key", as_index=False).mean(numeric_only=True)
        prot_summary_id["protein_abundance_rank_by_id"] = rank_pct(prot_summary_id["protein_abundance_value"])
    else:
        prot_summary_id = pd.DataFrame(columns=["gene_id_key","protein_abundance_value","protein_abundance_rank_by_id"])

    prot_summary_sym = prot[prot["gene_symbol_key"].astype(str).str.len() > 0][["gene_symbol_key","protein_abundance_value"]].groupby("gene_symbol_key", as_index=False).mean(numeric_only=True)
    prot_summary_sym["protein_abundance_rank_by_symbol"] = rank_pct(prot_summary_sym["protein_abundance_value"])

    print("  mapped proteomics gene IDs:", prot_summary_id["gene_id_key"].nunique())
    print("  mapped proteomics symbols:", prot_summary_sym["gene_symbol_key"].nunique())

    # Merge by gene_id first, then symbol fallback
    df = lab.merge(
        prot_summary_id.rename(columns={"protein_abundance_value":"protein_abundance_value_by_id"}),
        on="gene_id_key",
        how="left"
    )
    df = df.merge(
        prot_summary_sym.rename(columns={"protein_abundance_value":"protein_abundance_value_by_symbol"}),
        on="gene_symbol_key",
        how="left"
    )

    df["protein_abundance_value"] = df["protein_abundance_value_by_id"].fillna(df["protein_abundance_value_by_symbol"])
    df["protein_abundance_rank"] = df["protein_abundance_rank_by_id"].fillna(df["protein_abundance_rank_by_symbol"])
    df["proteomics_mapping_mode"] = np.where(df["protein_abundance_value_by_id"].notna(), "gene_id", np.where(df["protein_abundance_value_by_symbol"].notna(), "gene_symbol", "unmapped"))
    df["has_proteomics_label"] = df["protein_abundance_value"].notna()

    pc = 1e-6
    if "rna_day3_cpm_mean" in df.columns and "rna_day6_cpm_mean" in df.columns:
        avg_rna = ((df["rna_day3_cpm_mean"] + df["rna_day6_cpm_mean"]) / 2).fillna(0).values
    else:
        avg_rna = np.zeros(len(df))
    df["protein_residual_rank"] = residual_rank(np.log2(avg_rna + pc), df["protein_abundance_value"].values)
    df.loc[df["protein_abundance_value"].isna(), "protein_residual_rank"] = np.nan

    w = cfg["multi_omics_score_weights"]
    fallback = df["robust_public_te_rank"].fillna(0)
    pab = df["protein_abundance_rank"].fillna(fallback)
    pres = df["protein_residual_rank"].fillna(fallback)

    df["multi_omics_utr_score"] = (
        w.get("robust_public_te_rank", 0.50) * df["robust_public_te_rank"].fillna(0)
        + w.get("protein_residual_rank", 0.20) * pres
        + w.get("protein_abundance_rank", 0.15) * pab
        + w.get("day_consensus_TE_rank", 0.10) * df["day_consensus_TE_rank"].fillna(0)
        + w.get("tss_confidence_score", 0.05) * df["tss_confidence_score"].fillna(0)
    )
    df["multi_omics_utr_rank"] = rank_pct(df["multi_omics_utr_score"])

    cf = cfg["candidate_filter"]
    seq_len_col = "utr5_length_final" if "utr5_length_final" in df.columns else "utr5_length_tss_corrected"
    train_ready = df[seq_len_col].between(cf["training_min_length"], cf["training_max_length"]) & df["multi_omics_utr_rank"].notna() & df["gc_content"].between(cf["gc_min"], cf["gc_max"]) & (df["uaug_count"] <= cf["uaug_max"])
    select_ready = df[seq_len_col].between(cf["selection_min_length"], cf["selection_max_length"]) & df["multi_omics_utr_rank"].notna() & df["gc_content"].between(cf["gc_min"], cf["gc_max"]) & (df["uaug_count"] <= cf["uaug_max"])
    df["multiomics_training_ready_40_200bp"] = train_ready
    df["multiomics_selection_ready_50_100bp"] = select_ready

    df.to_csv(OUT_LABEL, index=False)
    df[df["multiomics_training_ready_40_200bp"]].to_csv(OUT_READY_40_200, index=False)
    df[df["multiomics_selection_ready_50_100bp"]].to_csv(OUT_READY_50_100, index=False)

    lines = [
        "Proteomics / multi-omics mapping summary - NCBI gene-id aware",
        "="*90,
        f"proteomics_file: {prot_path}",
        f"proteomics_shape: {prot.shape}",
        f"abundance_columns: {ab_cols}",
        f"proteomics_gene_ids: {prot_summary_id['gene_id_key'].nunique()}",
        f"proteomics_symbols: {prot_summary_sym['gene_symbol_key'].nunique()}",
        f"UTR_rows: {len(df)}",
        f"UTR_rows_with_proteomics_label: {int(df['has_proteomics_label'].sum())}",
        f"UTR_rows_mapped_by_gene_id: {int((df['proteomics_mapping_mode']=='gene_id').sum())}",
        f"UTR_rows_mapped_by_gene_symbol: {int((df['proteomics_mapping_mode']=='gene_symbol').sum())}",
        f"UTR_rows_training_ready_40_200: {int(train_ready.sum())}",
        f"UTR_rows_selection_ready_50_100: {int(select_ready.sum())}",
        "",
        "[mapping mode counts]",
        df["proteomics_mapping_mode"].value_counts(dropna=False).to_string(),
        "",
        "[protein abundance rank]",
        df["protein_abundance_rank"].describe().to_string(),
        "",
        "[protein residual rank]",
        df["protein_residual_rank"].describe().to_string(),
        "",
        "[multi omics UTR rank]",
        df["multi_omics_utr_rank"].describe().to_string(),
    ]
    OUT_PROT_SUMMARY.write_text("\n".join(lines), encoding="utf-8")
    print("[SAVED]", OUT_LABEL)
    print("[SAVED]", OUT_READY_40_200)
    print("[SAVED]", OUT_READY_50_100)
    print("[SAVED]", OUT_PROT_SUMMARY)
    print("\n".join(lines[:17]))

if __name__ == "__main__":
    main()
