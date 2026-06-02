from __future__ import annotations

from pathlib import Path
from collections import defaultdict, Counter
import argparse
import random
import re
import numpy as np
import pandas as pd

BASE = Path.cwd()
SEQ = "utr5_sequence_tss_corrected"
DEFAULT_INPUTS = [
    BASE / "04_te_labeling/tables/tss_corrected_5utr_multiomics_labels.csv",
    BASE / "04_te_labeling/tables/tss_corrected_5utr_robust_public_te_labels.csv",
]
OUT_TABLE = BASE / "04_te_labeling/tables/tss_corrected_5utr_with_seq_clusters.csv"
OUT_CLUSTER = BASE / "04_te_labeling/qc/jaccard_sequence_cluster_summary.csv"
OUT_REPORT = BASE / "04_te_labeling/qc/jaccard_sequence_cluster_report.txt"
for p in [OUT_TABLE.parent, OUT_CLUSTER.parent, OUT_REPORT.parent]:
    p.mkdir(parents=True, exist_ok=True)

BASES = "ACGT"
BASE_TO_INT = {b: i for i, b in enumerate(BASES)}


def clean_seq(x: object) -> str:
    if pd.isna(x):
        return ""
    return re.sub(r"[^ACGTN]", "", str(x).upper().replace("U", "T"))


def gc_content(seq: str) -> float:
    seq = clean_seq(seq)
    return (seq.count("G") + seq.count("C")) / len(seq) if seq else np.nan


def choose_input(path_arg: str | None = None) -> Path:
    if path_arg:
        p = Path(path_arg)
        if not p.exists():
            raise SystemExit(f"Input not found: {p}")
        return p
    for p in DEFAULT_INPUTS:
        if p.exists():
            return p
    raise SystemExit("No label table found. Run publicTE/multiomics pipeline first.")


def kmer_id(kmer: str) -> int | None:
    v = 0
    for b in kmer:
        if b not in BASE_TO_INT:
            return None
        v = v * 4 + BASE_TO_INT[b]
    return v


def kmer_set(seq: str, k: int) -> frozenset[int]:
    seq = clean_seq(seq)
    if len(seq) < k:
        return frozenset()
    out = set()
    for i in range(len(seq) - k + 1):
        kid = kmer_id(seq[i:i+k])
        if kid is not None:
            out.add(kid)
    return frozenset(out)


def jaccard(a: frozenset[int], b: frozenset[int]) -> float:
    if not a and not b:
        return 1.0
    if not a or not b:
        return 0.0
    inter = len(a & b)
    union = len(a | b)
    return inter / union if union else 0.0


def containment(a: frozenset[int], b: frozenset[int]) -> float:
    if not a and not b:
        return 1.0
    if not a or not b:
        return 0.0
    return len(a & b) / min(len(a), len(b))


class UnionFind:
    def __init__(self, n: int):
        self.parent = list(range(n))
        self.rank = [0] * n

    def find(self, x: int) -> int:
        while self.parent[x] != x:
            self.parent[x] = self.parent[self.parent[x]]
            x = self.parent[x]
        return x

    def union(self, a: int, b: int) -> None:
        ra, rb = self.find(a), self.find(b)
        if ra == rb:
            return
        if self.rank[ra] < self.rank[rb]:
            ra, rb = rb, ra
        self.parent[rb] = ra
        if self.rank[ra] == self.rank[rb]:
            self.rank[ra] += 1


def minhash_signature(kset: frozenset[int], seeds: list[tuple[int, int]], prime: int = 1000003) -> tuple[int, ...]:
    if not kset:
        return tuple([prime] * len(seeds))
    vals = []
    arr = list(kset)
    for a, b in seeds:
        vals.append(min(((a * x + b) % prime) for x in arr))
    return tuple(vals)


def add_basic_columns(df: pd.DataFrame) -> pd.DataFrame:
    df = df.copy()
    df[SEQ] = df[SEQ].map(clean_seq)
    if "utr5_length_final" in df.columns:
        df["length_for_cluster"] = pd.to_numeric(df["utr5_length_final"], errors="coerce")
    elif "utr5_length_tss_corrected" in df.columns:
        df["length_for_cluster"] = pd.to_numeric(df["utr5_length_tss_corrected"], errors="coerce")
    else:
        df["length_for_cluster"] = df[SEQ].str.len()
    if "gc_content" not in df.columns:
        df["gc_content"] = df[SEQ].map(gc_content)
    if "uaug_count" not in df.columns:
        df["uaug_count"] = df[SEQ].str.count("ATG")
    return df


def main():
    ap = argparse.ArgumentParser(description="Exact duplicate + k-mer Jaccard/containment sequence clustering for CHO 5'UTR candidates")
    ap.add_argument("--input", default=None)
    ap.add_argument("--k", type=int, default=6)
    ap.add_argument("--jaccard-threshold", type=float, default=0.85)
    ap.add_argument("--containment-threshold", type=float, default=0.90)
    ap.add_argument("--num-perm", type=int, default=80)
    ap.add_argument("--bands", type=int, default=20)
    ap.add_argument("--cluster-scope", choices=["all", "train_20_500", "selection_50_100"], default="all")
    ap.add_argument("--max-bucket-pairs", type=int, default=250000, help="safety cap for very large LSH buckets")
    args = ap.parse_args()

    path = choose_input(args.input)
    print("[LOAD]", path)
    df = add_basic_columns(pd.read_csv(path))
    n = len(df)

    if args.cluster_scope == "train_20_500":
        scope = df["length_for_cluster"].between(20, 500) & df["gc_content"].between(0.25, 0.80) & (df["uaug_count"] <= 3)
    elif args.cluster_scope == "selection_50_100":
        scope = df["length_for_cluster"].between(50, 100) & df["gc_content"].between(0.30, 0.75) & (df["uaug_count"] <= 1)
    else:
        scope = pd.Series(True, index=df.index)

    scope &= df[SEQ].str.len().ge(args.k)
    scope &= ~df[SEQ].str.contains("N", regex=False)
    scope_idx = np.where(scope.values)[0]
    print(f"[SCOPE] {args.cluster_scope}: {len(scope_idx):,}/{n:,} rows")

    uf = UnionFind(n)

    # exact duplicate union over all rows with valid sequence
    seq_to_indices: dict[str, list[int]] = defaultdict(list)
    for i, s in enumerate(df[SEQ].tolist()):
        if s:
            seq_to_indices[s].append(i)
    exact_groups = 0
    for idxs in seq_to_indices.values():
        if len(idxs) > 1:
            exact_groups += 1
            first = idxs[0]
            for j in idxs[1:]:
                uf.union(first, j)
    print(f"[EXACT] duplicate sequence groups: {exact_groups:,}")

    # LSH only within scope
    ksets = {}
    for i in scope_idx:
        ksets[i] = kmer_set(df.at[i, SEQ], args.k)

    rnd = random.Random(42)
    seeds = [(rnd.randint(1, 999983), rnd.randint(0, 999983)) for _ in range(args.num_perm)]
    if args.num_perm % args.bands != 0:
        raise SystemExit("num-perm must be divisible by bands")
    rows_per_band = args.num_perm // args.bands

    buckets: dict[tuple[int, tuple[int, ...]], list[int]] = defaultdict(list)
    for i in scope_idx:
        sig = minhash_signature(ksets[i], seeds)
        for b in range(args.bands):
            key = (b, sig[b*rows_per_band:(b+1)*rows_per_band])
            buckets[key].append(i)

    candidate_pairs = set()
    skipped_large = 0
    for members in buckets.values():
        m = len(members)
        if m < 2:
            continue
        possible = m * (m - 1) // 2
        if possible > args.max_bucket_pairs:
            skipped_large += 1
            # sample deterministic subset of pairs to avoid blowup
            members = sorted(members)[:int((2 * args.max_bucket_pairs) ** 0.5)]
            m = len(members)
        for a_pos in range(m):
            a = members[a_pos]
            for b_pos in range(a_pos + 1, m):
                b = members[b_pos]
                candidate_pairs.add((a, b) if a < b else (b, a))

    print(f"[LSH] candidate pairs: {len(candidate_pairs):,}; skipped large buckets: {skipped_large}")

    linked = 0
    max_j = 0.0
    max_c = 0.0
    for a, b in candidate_pairs:
        ja = jaccard(ksets[a], ksets[b])
        co = containment(ksets[a], ksets[b])
        max_j = max(max_j, ja)
        max_c = max(max_c, co)
        if ja >= args.jaccard_threshold or co >= args.containment_threshold:
            uf.union(a, b)
            linked += 1
    print(f"[LINK] linked pairs: {linked:,}; max_jaccard={max_j:.3f}; max_containment={max_c:.3f}")

    roots = [uf.find(i) for i in range(n)]
    root_to_cluster = {}
    cluster_ids = []
    for r in roots:
        if r not in root_to_cluster:
            root_to_cluster[r] = f"SC{len(root_to_cluster)+1:06d}"
        cluster_ids.append(root_to_cluster[r])
    df["seq_cluster_id"] = cluster_ids

    sizes = df.groupby("seq_cluster_id").size().rename("seq_cluster_size")
    df = df.merge(sizes, left_on="seq_cluster_id", right_index=True, how="left")
    df["exact_duplicate_count"] = df[SEQ].map(df.groupby(SEQ).size())

    # Representative: highest robust_public_te_rank if present, else first.
    sort_cols = []
    asc = []
    for c in ["robust_public_te_rank", "multi_omics_utr_rank", "protein_residual_rank", "protein_abundance_rank"]:
        if c in df.columns:
            sort_cols.append(c)
            asc.append(False)
    sort_cols.append("length_for_cluster")
    asc.append(True)
    rep_idx = df.sort_values(sort_cols, ascending=asc).groupby("seq_cluster_id", as_index=False).head(1).set_index("seq_cluster_id").index
    reps = df.sort_values(sort_cols, ascending=asc).groupby("seq_cluster_id", as_index=False).head(1)[["seq_cluster_id", SEQ]]
    reps = reps.rename(columns={SEQ: "seq_cluster_representative_sequence"})
    df = df.merge(reps, on="seq_cluster_id", how="left")

    cluster_rows = []
    rank_cols = [c for c in ["robust_public_te_rank", "multi_omics_utr_rank", "protein_abundance_rank", "protein_residual_rank"] if c in df.columns]
    for cid, g in df.groupby("seq_cluster_id"):
        row = {
            "seq_cluster_id": cid,
            "size": len(g),
            "unique_genes": g["gene_name"].nunique() if "gene_name" in g.columns else np.nan,
            "length_min": g["length_for_cluster"].min(),
            "length_max": g["length_for_cluster"].max(),
            "representative_sequence": g["seq_cluster_representative_sequence"].iloc[0],
        }
        for c in rank_cols:
            row[f"{c}_mean"] = pd.to_numeric(g[c], errors="coerce").mean()
            row[f"{c}_max"] = pd.to_numeric(g[c], errors="coerce").max()
        cluster_rows.append(row)
    clus = pd.DataFrame(cluster_rows).sort_values("size", ascending=False)

    df.to_csv(OUT_TABLE, index=False)
    clus.to_csv(OUT_CLUSTER, index=False)

    report = [
        "Jaccard / containment sequence clustering report",
        "=" * 90,
        f"input: {path}",
        f"rows_total: {n}",
        f"scope: {args.cluster_scope}",
        f"scope_rows: {len(scope_idx)}",
        f"k: {args.k}",
        f"jaccard_threshold: {args.jaccard_threshold}",
        f"containment_threshold: {args.containment_threshold}",
        f"num_perm: {args.num_perm}",
        f"bands: {args.bands}",
        f"exact_duplicate_sequence_groups: {exact_groups}",
        f"candidate_pairs_lsh: {len(candidate_pairs)}",
        f"linked_pairs: {linked}",
        f"clusters_total: {df['seq_cluster_id'].nunique()}",
        f"clusters_size_gt1: {(sizes > 1).sum()}",
        "",
        "[Cluster size summary]",
        sizes.describe().to_string(),
        "",
        "[Top clusters]",
        clus.head(30).to_string(index=False),
        "",
        f"Saved table: {OUT_TABLE}",
        f"Saved cluster summary: {OUT_CLUSTER}",
    ]
    OUT_REPORT.write_text("\n".join(report), encoding="utf-8")
    print("[SAVED]", OUT_TABLE)
    print("[SAVED]", OUT_CLUSTER)
    print("[SAVED]", OUT_REPORT)


if __name__ == "__main__":
    main()
