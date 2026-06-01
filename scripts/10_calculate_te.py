import argparse
from pathlib import Path
import pandas as pd
import numpy as np


def safe_log2(x):
    return np.log2(x + 1e-9)


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--config", default="config/config.yaml")
    parser.add_argument("--rna", default="data/processed/rna_abundance.csv")
    parser.add_argument("--ribo", default="data/processed/ribo_abundance.csv")
    parser.add_argument("--protein", default="data/processed/protein_abundance.csv")
    parser.add_argument("--output", default="data/processed/te_metrics.csv")
    args = parser.parse_args()

    rna = pd.read_csv(args.rna)
    ribo = pd.read_csv(args.ribo)
    protein = pd.read_csv(args.protein)

    required_rna = {"gene_id", "rna_abundance"}
    required_ribo = {"gene_id", "ribo_abundance"}
    required_protein = {"gene_id", "protein_abundance"}

    if not required_rna.issubset(rna.columns):
        raise ValueError(f"RNA file must contain columns: {required_rna}")
    if not required_ribo.issubset(ribo.columns):
        raise ValueError(f"Ribo file must contain columns: {required_ribo}")
    if not required_protein.issubset(protein.columns):
        raise ValueError(f"Protein file must contain columns: {required_protein}")

    df = rna.merge(ribo, on="gene_id", how="inner")
    df = df.merge(protein, on="gene_id", how="inner")

    df["log2_rna"] = safe_log2(df["rna_abundance"])
    df["log2_ribo"] = safe_log2(df["ribo_abundance"])
    df["log2_protein"] = safe_log2(df["protein_abundance"])

    df["ribo_te"] = df["log2_ribo"] - df["log2_rna"]
    df["protein_te"] = df["log2_protein"] - df["log2_rna"]

    # Protein residual rank:
    # protein이 높은 이유가 단순히 RNA가 많아서인지 제거하기 위해
    # log2(protein) ~ log2(RNA) 선형 보정 후 residual 계산
    x = df["log2_rna"].values
    y = df["log2_protein"].values

    slope, intercept = np.polyfit(x, y, 1)
    predicted = slope * x + intercept

    df["protein_predicted_from_rna"] = predicted
    df["protein_residual"] = y - predicted

    df["ribo_te_rank"] = df["ribo_te"].rank(method="average", ascending=False)
    df["protein_te_rank"] = df["protein_te"].rank(method="average", ascending=False)
    df["protein_residual_rank"] = df["protein_residual"].rank(method="average", ascending=False)

    out = Path(args.output)
    out.parent.mkdir(parents=True, exist_ok=True)
    df.to_csv(out, index=False)

    print(f"[OK] TE metrics saved to: {out}")
    print(f"[INFO] Number of matched genes: {len(df)}")


if __name__ == "__main__":
    main()