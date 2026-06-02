import argparse
from pathlib import Path
import pandas as pd


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--config", default="config/config.yaml")
    parser.add_argument("--input", default="data/processed/riboseq_counts_mapped.csv")
    parser.add_argument("--totals", default="data/processed/riboseq_total_mapped_reads.csv")
    parser.add_argument("--output", default="data/processed/ribo_abundance.csv")
    args = parser.parse_args()

    inp = Path(args.input)
    totals_path = Path(args.totals)

    if not inp.exists():
        raise FileNotFoundError(f"Missing mapped Ribo-seq counts: {inp}")
    if not totals_path.exists():
        raise FileNotFoundError(f"Missing Ribo-seq total mapped reads: {totals_path}")

    df = pd.read_csv(inp)
    totals = pd.read_csv(totals_path).iloc[0].to_dict()

    sample_cols = [c for c in df.columns if c.startswith("s")]

    for c in sample_cols:
        total = float(totals[c])
        df[f"{c}_cpm"] = df[c] / total * 1_000_000

    cpm_cols = [f"{c}_cpm" for c in sample_cols]

    df["ribo_abundance"] = df[cpm_cols].mean(axis=1)
    df["ribo_abundance_median"] = df[cpm_cols].median(axis=1)
    df["ribo_raw_count_mean"] = df[sample_cols].mean(axis=1)

    keep = [
        "gene_id",
        "gene_symbol",
        "transcript_id",
        "protein_id",
        "cds_nucleotide_length",
        "ribo_abundance",
        "ribo_abundance_median",
        "ribo_raw_count_mean",
    ] + sample_cols + cpm_cols

    out = Path(args.output)
    out.parent.mkdir(parents=True, exist_ok=True)
    df[keep].to_csv(out, index=False)

    print(f"[OK] Ribo abundance saved to: {out}")
    print(f"[INFO] Ribo abundance rows: {len(df)}")
    print(f"[INFO] CPM sample columns: {cpm_cols}")


if __name__ == "__main__":
    main()