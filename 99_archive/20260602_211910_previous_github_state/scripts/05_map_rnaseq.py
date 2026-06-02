import argparse
from pathlib import Path

import pandas as pd


def read_public_count_table(path: Path) -> tuple[pd.DataFrame, dict]:
    df = pd.read_csv(path, sep="\t", compression="infer")
    total = df[df["geneid"].astype(str) == "TotalMappedReads"]
    if total.empty:
        raise ValueError("TotalMappedReads row not found")
    sample_cols = [c for c in df.columns if c.startswith("s")]
    totals = {c: float(total.iloc[0][c]) for c in sample_cols}
    df = df[df["geneid"].astype(str) != "TotalMappedReads"].copy()
    for c in sample_cols:
        df[c] = pd.to_numeric(df[c], errors="coerce").fillna(0)
    return df, totals


def main():
    parser = argparse.ArgumentParser(description="Load GSE79512 RNA-seq raw count table and standardize identifiers.")
    parser.add_argument("--config", default="config/config.yaml")
    parser.add_argument("--input", default="data/raw/rna_ribo/GSE79512_RNASeq_rawCount.txt.gz")
    parser.add_argument("--output", default="data/processed/rnaseq_counts_mapped.csv")
    parser.add_argument("--total-output", default="data/processed/rnaseq_total_mapped_reads.csv")
    args = parser.parse_args()

    inp = Path(args.input)
    if not inp.exists():
        raise FileNotFoundError(f"Missing RNA-seq count file: {inp}")

    df, totals = read_public_count_table(inp)
    df = df.rename(columns={
        "geneid": "gene_id",
        "genesymbol": "gene_symbol",
        "rnaid": "transcript_id",
        "prid": "protein_id",
    })

    sample_cols = [c for c in df.columns if c.startswith("s")]
    keep_cols = ["gene_id", "gene_symbol", "transcript_id", "protein_id", "rna_length"] + sample_cols
    df = df[keep_cols]

    out = Path(args.output)
    total_out = Path(args.total_output)
    out.parent.mkdir(parents=True, exist_ok=True)
    total_out.parent.mkdir(parents=True, exist_ok=True)
    df.to_csv(out, index=False)
    pd.DataFrame([totals]).to_csv(total_out, index=False)

    print(f"[OK] RNA-seq mapped counts saved to: {out}")
    print(f"[OK] RNA-seq total mapped reads saved to: {total_out}")
    print(f"[INFO] RNA-seq transcripts: {len(df)}")
    print(f"[INFO] Sample columns: {sample_cols}")


if __name__ == "__main__":
    main()
