import argparse
from pathlib import Path
import pandas as pd
import numpy as np
import pickle


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


def make_features(df):
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
    return df


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--config", default="config/config.yaml")
    parser.add_argument("--input", default="data/candidates/candidate_utr_sequences.csv")
    parser.add_argument("--model", default="results/14_train_rank_model/model.pkl")
    parser.add_argument("--output", default="results/15_score_candidate_library/candidate_scores.csv")
    args = parser.parse_args()

    df = pd.read_csv(args.input)

    required = {"candidate_id", "utr_sequence"}
    if not required.issubset(df.columns):
        raise ValueError(f"Candidate file must contain columns: {required}")

    df = make_features(df)

    feature_cols = [
        "utr_length",
        "gc_content",
        "at_content",
        "num_a",
        "num_c",
        "num_g",
        "num_t",
        "num_uaug",
        "has_uaug",
        "kozak_proxy_score",
    ]

    model_path = Path(args.model)

    if model_path.exists():
        with open(model_path, "rb") as f:
            model = pickle.load(f)
        df["predicted_te_rank"] = model.predict(df[feature_cols])
        print(f"[INFO] Model loaded: {model_path}")
    else:
        # fallback heuristic until final model object is available
        df["predicted_te_rank"] = (
            50
            + 20 * df["gc_content"].fillna(0)
            - 5 * df["has_uaug"].fillna(0)
            + 2 * df["kozak_proxy_score"].fillna(0)
        )
        print("[WARN] Model file not found. Used temporary heuristic scoring.")

    df["candidate_rank"] = df["predicted_te_rank"].rank(method="first", ascending=False).astype(int)

    out = Path(args.output)
    out.parent.mkdir(parents=True, exist_ok=True)
    df.sort_values("candidate_rank").to_csv(out, index=False)

    print(f"[OK] Candidate scores saved to: {out}")
    print(f"[INFO] Number of candidates: {len(df)}")


if __name__ == "__main__":
    main()