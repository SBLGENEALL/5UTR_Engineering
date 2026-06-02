import argparse
import itertools
import re
import shutil
import subprocess
from pathlib import Path

import numpy as np
import pandas as pd


DNA_ALPHABET = "ACGT"


def clean_seq(seq):
    return re.sub(r"[^ACGT]", "N", str(seq).upper().replace("U", "T"))


def gc_content(seq):
    seq = clean_seq(seq)
    if not seq:
        return np.nan
    valid = [b for b in seq if b in DNA_ALPHABET]
    if not valid:
        return np.nan
    return (valid.count("G") + valid.count("C")) / len(valid)


def at_content(seq):
    seq = clean_seq(seq)
    if not seq:
        return np.nan
    valid = [b for b in seq if b in DNA_ALPHABET]
    if not valid:
        return np.nan
    return (valid.count("A") + valid.count("T")) / len(valid)


def count_motif(seq, motif):
    return clean_seq(seq).count(motif.upper().replace("U", "T"))


def has_uaug(seq):
    return int("ATG" in clean_seq(seq))


def kozak_score(seq):
    seq = clean_seq(seq)
    idx = seq.find("ATG")
    if idx < 0:
        return 0
    score = 0
    if idx - 3 >= 0 and seq[idx - 3] in ["A", "G"]:
        score += 1
    if idx + 3 < len(seq) and seq[idx + 3] == "G":
        score += 1
    return score


def make_kmers(kmax):
    kmers = []
    for k in range(1, kmax + 1):
        kmers.extend("".join(p) for p in itertools.product(DNA_ALPHABET, repeat=k))
    return kmers


def add_kmer_features(df, seq_col, kmax):
    if kmax <= 0:
        return df
    kmers = make_kmers(kmax)
    seqs = df[seq_col].map(clean_seq)
    for kmer in kmers:
        k = len(kmer)
        col = f"kmer_{kmer}"
        vals = []
        for seq in seqs:
            denom = max(len(seq) - k + 1, 1)
            vals.append(seq.count(kmer) / denom)
        df[col] = vals
    return df


def run_rnafold_one(seq, rnafold_bin="RNAfold"):
    seq = clean_seq(seq)
    if not seq or "N" in seq:
        return np.nan, ""
    proc = subprocess.run(
        [rnafold_bin, "--noPS"],
        input=seq + "\n",
        text=True,
        capture_output=True,
        check=False,
    )
    if proc.returncode != 0:
        return np.nan, ""
    lines = [x.strip() for x in proc.stdout.splitlines() if x.strip()]
    if len(lines) < 2:
        return np.nan, ""
    struct_line = lines[1]
    m = re.search(r"\(([-+]?\d+(?:\.\d+)?)\)", struct_line)
    mfe = float(m.group(1)) if m else np.nan
    structure = struct_line.split()[0]
    return mfe, structure


def add_rnafold_features(df, seq_col, rnafold_bin="RNAfold", head_tail=30):
    if shutil.which(rnafold_bin) is None:
        raise FileNotFoundError(f"RNAfold executable not found: {rnafold_bin}")

    full_mfe = []
    full_structure = []
    head_mfe = []
    tail_mfe = []

    seqs = df[seq_col].map(clean_seq).tolist()
    total = len(seqs)
    for i, seq in enumerate(seqs, start=1):
        if i == 1 or i % 500 == 0 or i == total:
            print(f"[RNAfold] {i}/{total}")
        mfe, struct = run_rnafold_one(seq, rnafold_bin=rnafold_bin)
        hmfe, _ = run_rnafold_one(seq[:head_tail], rnafold_bin=rnafold_bin)
        tmfe, _ = run_rnafold_one(seq[-head_tail:], rnafold_bin=rnafold_bin)
        full_mfe.append(mfe)
        full_structure.append(struct)
        head_mfe.append(hmfe)
        tail_mfe.append(tmfe)

    df["mfe"] = full_mfe
    df["mfe_per_nt"] = df["mfe"] / df["utr_length"].replace(0, np.nan)
    df[f"head{head_tail}_mfe"] = head_mfe
    df[f"head{head_tail}_mfe_per_nt"] = df[f"head{head_tail}_mfe"] / np.minimum(df["utr_length"], head_tail).replace(0, np.nan)
    df[f"tail{head_tail}_mfe"] = tail_mfe
    df[f"tail{head_tail}_mfe_per_nt"] = df[f"tail{head_tail}_mfe"] / np.minimum(df["utr_length"], head_tail).replace(0, np.nan)
    df["rnafold_structure"] = full_structure
    df["rnafold_paired_fraction"] = [s.count("(") * 2 / len(s) if s else np.nan for s in full_structure]
    return df


def main():
    parser = argparse.ArgumentParser(description="Extract final 5UTR features for model training.")
    parser.add_argument("--config", default="config/config.yaml")
    parser.add_argument("--input", default="data/processed/te_rank_labels.csv")
    parser.add_argument("--utr", default="data/processed/utr_sequences.csv")
    parser.add_argument("--output", default="data/processed/feature_matrix.csv")
    parser.add_argument("--min-length", type=int, default=40)
    parser.add_argument("--max-length", type=int, default=200)
    parser.add_argument("--kmax", type=int, default=4)
    parser.add_argument("--use-rnafold", action="store_true")
    parser.add_argument("--rnafold-bin", default="RNAfold")
    parser.add_argument("--head-tail", type=int, default=30)
    args = parser.parse_args()

    labels = pd.read_csv(args.input)
    utr = pd.read_csv(args.utr)

    if "gene_id" not in labels.columns:
        raise ValueError("Label file must contain gene_id")
    if "gene_id" not in utr.columns or "utr_sequence" not in utr.columns:
        raise ValueError("UTR file must contain gene_id and utr_sequence")

    labels["gene_id"] = labels["gene_id"].astype(str)
    utr["gene_id"] = utr["gene_id"].astype(str)
    utr["utr_sequence"] = utr["utr_sequence"].map(clean_seq)

    df = labels.merge(utr, on="gene_id", how="inner", suffixes=("", "_utr"))
    df["utr_length"] = df["utr_sequence"].astype(str).str.len()
    df = df[(df["utr_status"] == "ok") & (df["utr_length"] >= args.min_length) & (df["utr_length"] <= args.max_length)].copy()
    df = df.reset_index(drop=True)

    df["utr5_sequence_tss_corrected"] = df["utr_sequence"]
    df["length_for_cluster"] = df["utr_length"]
    df["gc_content"] = df["utr_sequence"].apply(gc_content)
    df["at_content"] = df["utr_sequence"].apply(at_content)
    df["num_a"] = df["utr_sequence"].str.count("A")
    df["num_c"] = df["utr_sequence"].str.count("C")
    df["num_g"] = df["utr_sequence"].str.count("G")
    df["num_t"] = df["utr_sequence"].str.count("T")
    df["num_n"] = df["utr_sequence"].str.count("N")
    df["num_uaug"] = df["utr_sequence"].apply(lambda x: count_motif(x, "ATG"))
    df["uaug_count"] = df["num_uaug"]
    df["has_uaug"] = df["utr_sequence"].apply(has_uaug)
    df["kozak_proxy_score"] = df["utr_sequence"].apply(kozak_score)

    df = add_kmer_features(df, "utr_sequence", args.kmax)

    if args.use_rnafold:
        df = add_rnafold_features(df, "utr_sequence", rnafold_bin=args.rnafold_bin, head_tail=args.head_tail)
    else:
        df["mfe"] = np.nan
        df["mfe_per_nt"] = np.nan
        df[f"head{args.head_tail}_mfe"] = np.nan
        df[f"head{args.head_tail}_mfe_per_nt"] = np.nan
        df[f"tail{args.head_tail}_mfe"] = np.nan
        df[f"tail{args.head_tail}_mfe_per_nt"] = np.nan
        df["rnafold_structure"] = ""
        df["rnafold_paired_fraction"] = np.nan

    out = Path(args.output)
    out.parent.mkdir(parents=True, exist_ok=True)
    df.to_csv(out, index=False)

    print(f"[OK] Feature matrix saved to: {out}")
    print(f"[INFO] Length window: {args.min_length}-{args.max_length}")
    print(f"[INFO] Rows: {len(df)}")
    print(f"[INFO] Columns: {df.shape[1]}")
    print(f"[INFO] kmax: {args.kmax}")
    print(f"[INFO] RNAfold: {args.use_rnafold}")


if __name__ == "__main__":
    main()
