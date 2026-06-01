import argparse
from pathlib import Path
import pandas as pd


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--config", default="config/config.yaml")
    parser.add_argument(
        "--input",
        default="data/raw/proteomics/Heffner_2020_CHO_hamster_proteomics_ncbi_mapped_for_5utr.csv",
    )
    parser.add_argument("--output", default="data/processed/protein_abundance.csv")
    args = parser.parse_args()

    inp = Path(args.input)
    if not inp.exists():
        raise FileNotFoundError(f"Missing proteomics file: {inp}")

    df = pd.read_csv(inp)

    required = {"gene_id_key", "gene_symbol", "abundance_proxy"}
    if not required.issubset(df.columns):
        raise ValueError(f"Input must contain columns: {required}")

    out_df = df.copy()
    out_df = out_df.rename(columns={
        "gene_id_key": "gene_id",
        "abundance_proxy": "protein_abundance",
    })

    out_df["gene_id"] = out_df["gene_id"].astype(str)
    out_df["protein_abundance"] = pd.to_numeric(
        out_df["protein_abundance"], errors="coerce"
    )

    out_df = out_df.dropna(subset=["gene_id", "protein_abundance"])

    out_df = (
        out_df
        .groupby(["gene_id", "gene_symbol"], as_index=False)
        .agg(
            protein_abundance=("protein_abundance", "mean"),
            protein_entry_count=("protein_abundance", "size"),
        )
    )

    out = Path(args.output)
    out.parent.mkdir(parents=True, exist_ok=True)
    out_df.to_csv(out, index=False)

    print(f"[OK] Protein abundance saved to: {out}")
    print(f"[INFO] Protein abundance rows: {len(out_df)}")


if __name__ == "__main__":
    main()