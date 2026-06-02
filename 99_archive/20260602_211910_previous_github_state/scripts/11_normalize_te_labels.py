import argparse
from pathlib import Path
import pandas as pd


def add_percentile_rank(df, value_col, out_col):
    """
    Higher value = better rank.
    Output range: 0 to 100.
    """
    df[out_col] = df[value_col].rank(pct=True, ascending=True) * 100
    return df


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--config", default="config/config.yaml")
    parser.add_argument("--input", default="data/processed/te_metrics.csv")
    parser.add_argument("--output", default="data/processed/te_rank_labels.csv")
    parser.add_argument(
        "--label-source",
        default="protein_residual",
        choices=["ribo_te", "protein_te", "protein_residual"],
        help="Metric used as the main ML label",
    )
    args = parser.parse_args()

    df = pd.read_csv(args.input)

    required = {"gene_id", args.label_source}
    if not required.issubset(df.columns):
        raise ValueError(f"Input must contain columns: {required}")

    df = add_percentile_rank(
        df,
        value_col=args.label_source,
        out_col="te_percentile_rank",
    )

    df["te_rank_label"] = df["te_percentile_rank"]

    df["te_class_20_500"] = "middle"
    df.loc[df["te_percentile_rank"] >= 80, "te_class_20_500"] = "top20"
    df.loc[df["te_percentile_rank"] <= 20, "te_class_20_500"] = "bottom20"

    df["te_class_50_100"] = "middle"
    df.loc[df["te_percentile_rank"] >= 50, "te_class_50_100"] = "top50"
    df.loc[df["te_percentile_rank"] <= 50, "te_class_50_100"] = "bottom50"

    out = Path(args.output)
    out.parent.mkdir(parents=True, exist_ok=True)
    df.to_csv(out, index=False)

    print(f"[OK] TE rank labels saved to: {out}")
    print(f"[INFO] Label source: {args.label_source}")
    print(f"[INFO] Number of rows: {len(df)}")


if __name__ == "__main__":
    main()