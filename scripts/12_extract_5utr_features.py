import argparse
from pathlib import Path
import pandas as pd
import numpy as np


def gc_content(seq):
    seq = str(seq).upper()
    if len(seq) == 0:
        return np.nan
    return (seq.count("G") + seq.count("C")) / len(seq)


def at_content(seq):
    seq = str(seq).upper()
    if len(seq) == 0:
        return np.nan
    return (seq.count("A") + seq.count("T") + seq.count("U")) / len(seq)


def count_motif(seq, motif):
    return str(seq).upper().count(motif.upper())


def has_uaug(seq):
    seq = str(seq).upper().replace("U", "T")
    return int("ATG" in seq)


def kozak_score(seq):
    """
    Simple start-context proxy.
    If the sequence includes ATG, score bases around the first ATG.
    Higher score roughly means stronger Kozak-like context.
    """
    seq = str(seq).upper().replace("U", "T")
    idx = seq.find("ATG")
    if idx < 0:
        return 0

    score = 0
    if idx - 3 >= 0 and seq[idx - 3] in ["A", "G"]:
        score += 1
    if idx + 3 < len(seq) and seq[idx + 3] == "G":
        score += 1
    return score


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--config", default="config/config.yaml")
    parser.add_argument("--input", default="data/processed/te_rank_labels.csv")
    parser.add_argument("--utr", default="data/processed/utr_sequences.csv")
    parser.add_argument("--output", default="data/processed/feature_matrix.csv")
    args = parser.parse_args()

    labels = pd.read_csv(args.input)
    utr = pd.read_csv(args.utr)

    required_utr = {"gene_id", "utr_sequence"}
    if not required_utr.issubset(utr.columns):
        raise ValueError(f"UTR file must contain columns: {required_utr}")

    df = labels.merge(utr, on="gene_id", how="inner")

    df["utr_length"] = df["utr_sequence"].astype(str).str.len()
    df["gc_content"] = df["utr_sequence"].apply(gc_content)
    df["at_content"] = df["utr_sequence"].apply(at_content)
    df["num_a"] = df["utr_sequence"].astype(str).str.upper().str.count("A")
    df["num_c"] = df["utr_sequence"].astype(str).str.upper().str.count("C")
    df["num_g"] = df["utr_sequence"].astype(str).str.upper().str.count("G")
    df["num_t"] = df["utr_sequence"].astype(str).str.upper().str.count("T")
    df["num_uaug"] = df["utr_sequence"].apply(lambda x: count_motif(str(x).replace("U", "T"), "ATG"))
    df["has_uaug"] = df["utr_sequence"].apply(has_uaug)
    df["kozak_proxy_score"] = df["utr_sequence"].apply(kozak_score)

    out = Path(args.output)
    out.parent.mkdir(parents=True, exist_ok=True)
    df.to_csv(out, index=False)

    print(f"[OK] Feature matrix saved to: {out}")
    print(f"[INFO] Number of rows: {len(df)}")
    print(f"[INFO] Number of columns: {df.shape[1]}")


if __name__ == "__main__":
    main()