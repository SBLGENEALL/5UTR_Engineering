from __future__ import annotations

from pathlib import Path
from collections import Counter
from itertools import product
import argparse
import re
import math
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

from sklearn.model_selection import train_test_split
from sklearn.pipeline import Pipeline
from sklearn.impute import SimpleImputer
from sklearn.ensemble import ExtraTreesClassifier, RandomForestClassifier, HistGradientBoostingClassifier
from sklearn.metrics import roc_auc_score, average_precision_score

BASE = Path.cwd()
SEQ = "utr5_sequence_tss_corrected"
DEFAULT_INPUTS = [
    BASE / "04_te_labeling/tables/tss_corrected_5utr_with_seq_clusters.csv",
    BASE / "04_te_labeling/tables/tss_corrected_5utr_multiomics_labels.csv",
    BASE / "04_te_labeling/tables/tss_corrected_5utr_robust_public_te_labels.csv",
]
OUT = BASE / "06_modeling/tables/cluster_aware_classification_benchmark.csv"
DISJOINT_OUT = BASE / "06_modeling/tables/cluster_split_disjointness_check.csv"
PLOT = BASE / "06_modeling/plots/cluster_aware_classification_benchmark.png"
REPORT = BASE / "06_modeling/tables/cluster_aware_classification_summary.txt"
for p in [OUT.parent, PLOT.parent]:
    p.mkdir(parents=True, exist_ok=True)


def clean_seq(x):
    if pd.isna(x):
        return ""
    return re.sub(r"[^ACGTN]", "", str(x).upper().replace("U", "T"))


def gc_content(seq):
    seq = clean_seq(seq)
    return (seq.count("G") + seq.count("C")) / len(seq) if seq else np.nan


def longest_run(seq, base):
    best = cur = 0
    for b in clean_seq(seq):
        if b == base:
            cur += 1
            best = max(best, cur)
        else:
            cur = 0
    return best


def count_uorf(seq):
    seq = clean_seq(seq)
    stops = {"TAA", "TAG", "TGA"}
    n, longest = 0, 0
    for i in range(max(0, len(seq) - 2)):
        if seq[i:i+3] == "ATG":
            for j in range(i+3, len(seq)-2, 3):
                if seq[j:j+3] in stops:
                    n += 1
                    longest = max(longest, j+3-i)
                    break
    return n, longest


def kmer_features(seq, ks=(3,4,5,6), prefix=""):
    seq = clean_seq(seq)
    out = {}
    alphabet = "ACGT"
    for k in ks:
        all_k = ["".join(x) for x in product(alphabet, repeat=k)]
        cnt = Counter(seq[i:i+k] for i in range(max(0, len(seq)-k+1)) if set(seq[i:i+k]) <= set(alphabet))
        total = sum(cnt.values())
        for km in all_k:
            out[f"{prefix}{k}mer_{km}"] = cnt.get(km, 0) / total if total else 0.0
    return out


def make_features(df, kmax=5):
    rows = []
    ks = tuple(range(3, kmax+1))
    for seq in df[SEQ].map(clean_seq).tolist():
        L = len(seq)
        tail100 = seq[-100:] if L > 100 else seq
        head100 = seq[:100]
        uorf, longest = count_uorf(seq)
        r = {
            "length": L,
            "gc_content_feat": gc_content(seq),
            "tail100_gc": gc_content(tail100),
            "head100_gc": gc_content(head100),
            "uaug_count_feat": seq.count("ATG"),
            "tail100_uaug": tail100.count("ATG"),
            "head100_uaug": head100.count("ATG"),
            "cpg_count": seq.count("CG"),
            "first_aug_pos": seq.find("ATG"),
            "uorf_count": uorf,
            "longest_uorf_len": longest,
        }
        for b in "ACGT":
            r[f"{b.lower()}_fraction"] = seq.count(b) / L if L else 0.0
            r[f"poly{b}_max_run"] = longest_run(seq, b)
            r[f"tail100_{b.lower()}_fraction"] = tail100.count(b) / len(tail100) if len(tail100) else 0.0
        for w in [20, 30, 50, 100]:
            head = seq[:w]
            tail = seq[-w:] if L >= w else seq
            r[f"head{w}_gc"] = gc_content(head)
            r[f"tail{w}_gc"] = gc_content(tail)
            r[f"head{w}_uaug"] = head.count("ATG")
            r[f"tail{w}_uaug"] = tail.count("ATG")
        r.update(kmer_features(seq, ks=ks, prefix="full_"))
        r.update(kmer_features(tail100, ks=tuple(range(3, min(kmax,5)+1)), prefix="tail100_"))
        rows.append(r)
    return pd.DataFrame(rows)


class UnionFind:
    def __init__(self):
        self.parent = {}
        self.rank = {}
    def add(self, x):
        if x not in self.parent:
            self.parent[x] = x
            self.rank[x] = 0
    def find(self, x):
        self.add(x)
        while self.parent[x] != x:
            self.parent[x] = self.parent[self.parent[x]]
            x = self.parent[x]
        return x
    def union(self, a, b):
        self.add(a); self.add(b)
        ra, rb = self.find(a), self.find(b)
        if ra == rb:
            return
        if self.rank[ra] < self.rank[rb]:
            ra, rb = rb, ra
        self.parent[rb] = ra
        if self.rank[ra] == self.rank[rb]:
            self.rank[ra] += 1


def choose_input(path_arg=None):
    if path_arg:
        p = Path(path_arg)
        if p.exists(): return p
        raise SystemExit(f"Input not found: {p}")
    for p in DEFAULT_INPUTS:
        if p.exists(): return p
    raise SystemExit("No label table found")


def add_basic(df):
    df = df.copy()
    df[SEQ] = df[SEQ].map(clean_seq)
    if "length_for_cluster" in df.columns:
        df["length_for_model"] = pd.to_numeric(df["length_for_cluster"], errors="coerce")
    elif "utr5_length_final" in df.columns:
        df["length_for_model"] = pd.to_numeric(df["utr5_length_final"], errors="coerce")
    else:
        df["length_for_model"] = df[SEQ].str.len()
    if "gc_content" not in df.columns:
        df["gc_content"] = df[SEQ].map(gc_content)
    if "uaug_count" not in df.columns:
        df["uaug_count"] = df[SEQ].str.count("ATG")
    if "seq_cluster_id" not in df.columns:
        df["seq_cluster_id"] = "NOCLUSTER_" + df.index.astype(str)
    return df


def make_groups(df, split_mode):
    if split_mode == "random":
        return None
    if split_mode == "gene_split":
        return df["gene_name"].astype(str).values if "gene_name" in df.columns else df.index.astype(str).values
    if split_mode == "seq_cluster_split":
        return df["seq_cluster_id"].astype(str).values
    if split_mode == "gene_seq_cluster_split":
        uf = UnionFind()
        for idx, row in df.iterrows():
            rid = f"row:{idx}"
            uf.add(rid)
            if "gene_name" in df.columns and pd.notna(row.get("gene_name")):
                uf.union(rid, f"gene:{row['gene_name']}")
            uf.union(rid, f"cluster:{row['seq_cluster_id']}")
        return np.array([uf.find(f"row:{idx}") for idx in df.index])
    raise ValueError(split_mode)


def nonempty_values(df, col):
    if col not in df.columns:
        return set()
    vals = df[col].dropna().astype(str).str.strip()
    vals = vals[(vals != "") & (vals.str.lower() != "nan")]
    return set(vals)


def gene_key_column(df):
    for col in ["gene_name", "gene_id"]:
        if nonempty_values(df, col):
            return col
    return None


def overlap_examples(values, limit=10):
    values = sorted(values)
    return ";".join(values[:limit])


def disjointness_row(df, tr, te, target, split_mode):
    train = df.iloc[tr]
    test = df.iloc[te]
    gene_col = gene_key_column(df)
    if gene_col:
        train_genes = nonempty_values(train, gene_col)
        test_genes = nonempty_values(test, gene_col)
    else:
        train_genes = set()
        test_genes = set()

    train_clusters = nonempty_values(train, "seq_cluster_id")
    test_clusters = nonempty_values(test, "seq_cluster_id")
    gene_overlap = train_genes & test_genes
    cluster_overlap = train_clusters & test_clusters
    pass_gene = len(gene_overlap) == 0 if gene_col else False
    pass_cluster = len(cluster_overlap) == 0

    if split_mode == "gene_split":
        pass_required = pass_gene
    elif split_mode == "seq_cluster_split":
        pass_required = pass_cluster
    elif split_mode == "gene_seq_cluster_split":
        pass_required = pass_gene and pass_cluster
    else:
        pass_required = True

    return {
        "target": target,
        "split_mode": split_mode,
        "status": "ok",
        "n_total": len(df),
        "n_train": len(tr),
        "n_test": len(te),
        "gene_key_column": gene_col or "MISSING",
        "train_unique_genes": len(train_genes),
        "test_unique_genes": len(test_genes),
        "gene_overlap_count": len(gene_overlap),
        "gene_overlap_examples": overlap_examples(gene_overlap),
        "train_unique_seq_clusters": len(train_clusters),
        "test_unique_seq_clusters": len(test_clusters),
        "seq_cluster_overlap_count": len(cluster_overlap),
        "seq_cluster_overlap_examples": overlap_examples(cluster_overlap),
        "pass_gene_disjoint": pass_gene,
        "pass_seq_cluster_disjoint": pass_cluster,
        "pass_required_for_split": pass_required,
    }


def split_indices(df, y, split_mode):
    if split_mode == "random":
        return train_test_split(np.arange(len(df)), test_size=0.2, random_state=42, stratify=y)
    groups = make_groups(df, split_mode)
    uniq = np.array(sorted(pd.unique(groups)))
    # Try up to 100 random splits to retain both classes in test.
    for seed in range(42, 142):
        tr_g, te_g = train_test_split(uniq, test_size=0.2, random_state=seed)
        tr = np.where(np.isin(groups, tr_g))[0]
        te = np.where(np.isin(groups, te_g))[0]
        if len(np.unique(y[tr])) == 2 and len(np.unique(y[te])) == 2:
            return tr, te
    raise RuntimeError(f"Could not make valid split for {split_mode}")


def model_factory(name, n_estimators):
    if name == "ExtraTrees":
        return Pipeline([
            ("imputer", SimpleImputer(strategy="median")),
            ("model", ExtraTreesClassifier(n_estimators=n_estimators, max_features="sqrt", min_samples_leaf=3, class_weight="balanced", random_state=42, n_jobs=-1)),
        ])
    if name == "RandomForest":
        return Pipeline([
            ("imputer", SimpleImputer(strategy="median")),
            ("model", RandomForestClassifier(n_estimators=n_estimators, max_features="sqrt", min_samples_leaf=3, class_weight="balanced", random_state=42, n_jobs=-1)),
        ])
    if name == "HistGradientBoosting":
        return Pipeline([
            ("imputer", SimpleImputer(strategy="median")),
            ("model", HistGradientBoostingClassifier(max_iter=500, learning_rate=0.04, max_leaf_nodes=31, l2_regularization=0.02, random_state=42)),
        ])
    raise ValueError(name)


def main():
    ap = argparse.ArgumentParser(description="Cluster-aware high/low classification benchmark")
    ap.add_argument("--input", default=None)
    ap.add_argument("--length-min", type=int, default=20)
    ap.add_argument("--length-max", type=int, default=500)
    ap.add_argument("--kmax", type=int, default=5)
    ap.add_argument("--n-estimators", type=int, default=1000)
    ap.add_argument("--top-quantile", type=float, default=0.80)
    ap.add_argument("--bottom-quantile", type=float, default=0.30)
    args = ap.parse_args()

    path = choose_input(args.input)
    print("[LOAD]", path)
    df = add_basic(pd.read_csv(path))
    mask = df["length_for_model"].between(args.length_min, args.length_max) & df["gc_content"].between(0.25, 0.80) & (df["uaug_count"] <= 3)
    mask &= df[SEQ].str.len().gt(0) & ~df[SEQ].str.contains("N", regex=False)
    df = df[mask].drop_duplicates(subset=[SEQ]).reset_index(drop=True)
    print("[DATA]", df.shape)

    targets = [c for c in ["robust_public_te_rank", "multi_omics_utr_rank", "protein_abundance_rank", "protein_residual_rank"] if c in df.columns]
    split_modes = ["random", "gene_split", "seq_cluster_split", "gene_seq_cluster_split"]
    models = ["ExtraTrees", "RandomForest", "HistGradientBoosting"]

    X = make_features(df, kmax=args.kmax)
    rows = []
    disjoint_rows = []
    for target in targets:
        valid = pd.to_numeric(df[target], errors="coerce").notna()
        if valid.sum() < 500:
            print("[SKIP target few rows]", target, int(valid.sum()))
            continue
        sub = df[valid].copy().reset_index(drop=True)
        Xsub = X.loc[valid].reset_index(drop=True)
        yraw = pd.to_numeric(sub[target], errors="coerce").values
        hi = np.quantile(yraw, args.top_quantile)
        lo = np.quantile(yraw, args.bottom_quantile)
        keep = (yraw >= hi) | (yraw <= lo)
        sub = sub.loc[keep].copy().reset_index(drop=True)
        Xc = Xsub.loc[keep].reset_index(drop=True)
        y = (pd.to_numeric(sub[target], errors="coerce").values >= hi).astype(int)
        if len(sub) < 300 or len(np.unique(y)) < 2:
            print("[SKIP class few]", target, len(sub), np.bincount(y) if len(y) else [])
            continue
        print(f"[TARGET] {target}: n={len(sub)}, high={y.sum()}, low={(1-y).sum()}, hi_thr={hi:.3f}, lo_thr={lo:.3f}")

        for split_mode in split_modes:
            try:
                tr, te = split_indices(sub, y, split_mode)
            except Exception as e:
                rows.append({"target": target, "split_mode": split_mode, "model": "NA", "error": str(e)})
                if split_mode != "random":
                    disjoint_rows.append({
                        "target": target,
                        "split_mode": split_mode,
                        "status": f"split_failed: {e}",
                        "n_total": len(sub),
                        "n_train": np.nan,
                        "n_test": np.nan,
                        "gene_key_column": gene_key_column(sub) or "MISSING",
                        "train_unique_genes": np.nan,
                        "test_unique_genes": np.nan,
                        "gene_overlap_count": np.nan,
                        "gene_overlap_examples": "",
                        "train_unique_seq_clusters": np.nan,
                        "test_unique_seq_clusters": np.nan,
                        "seq_cluster_overlap_count": np.nan,
                        "seq_cluster_overlap_examples": "",
                        "pass_gene_disjoint": False,
                        "pass_seq_cluster_disjoint": False,
                        "pass_required_for_split": False,
                    })
                print("  [SPLIT FAIL]", split_mode, e)
                continue
            if split_mode != "random":
                disjoint_rows.append(disjointness_row(sub, tr, te, target, split_mode))
            for model_name in models:
                clf = model_factory(model_name, args.n_estimators)
                clf.fit(Xc.iloc[tr], y[tr])
                prob = clf.predict_proba(Xc.iloc[te])[:, 1]
                auc = roc_auc_score(y[te], prob)
                ap_score = average_precision_score(y[te], prob)
                rows.append({
                    "task": "classification",
                    "target": target,
                    "split_mode": split_mode,
                    "model": model_name,
                    "n_total": len(sub),
                    "n_train": len(tr),
                    "n_test": len(te),
                    "n_high_total": int(y.sum()),
                    "n_low_total": int((1-y).sum()),
                    "roc_auc": auc,
                    "average_precision": ap_score,
                    "selection_metric": 0.45 * auc + 0.55 * ap_score,
                    "top_quantile": args.top_quantile,
                    "bottom_quantile": args.bottom_quantile,
                    "length_min": args.length_min,
                    "length_max": args.length_max,
                    "kmax": args.kmax,
                    "n_estimators": args.n_estimators,
                })
                print(f"  {split_mode:24s} {model_name:22s} AUC={auc:.3f} AP={ap_score:.3f}")

    res = pd.DataFrame(rows)
    res.to_csv(OUT, index=False)
    disjoint = pd.DataFrame(disjoint_rows)
    disjoint.to_csv(DISJOINT_OUT, index=False)
    failed_disjoint = pd.DataFrame()
    if len(disjoint):
        failed_disjoint = disjoint[~disjoint["pass_required_for_split"].fillna(False).astype(bool)].copy()

    ok = res[res.get("roc_auc", pd.Series(dtype=float)).notna()].copy()
    if len(ok):
        sub = ok.sort_values("selection_metric", ascending=True).tail(30).copy()
        sub["label"] = sub["target"] + "\n" + sub["split_mode"] + "\n" + sub["model"]
        plt.figure(figsize=(12, max(6, 0.45 * len(sub))))
        plt.barh(sub["label"], sub["selection_metric"])
        plt.xlabel("0.45*ROC-AUC + 0.55*AveragePrecision")
        plt.title("Cluster-aware classification benchmark")
        plt.tight_layout()
        plt.savefig(PLOT, dpi=220)
        plt.close()

    lines = [
        "Cluster-aware classification benchmark summary",
        "=" * 90,
        f"input: {path}",
        f"rows after filter/dedup: {len(df)}",
        f"targets: {targets}",
        f"output: {OUT}",
        f"disjointness_check: {DISJOINT_OUT}",
        f"disjointness_required_failures: {len(failed_disjoint)}",
        f"plot: {PLOT}",
        "",
        "[Top rows by selection_metric]",
        ok.sort_values("selection_metric", ascending=False).head(40).to_string(index=False) if len(ok) else "No successful rows.",
    ]
    REPORT.write_text("\n".join(lines), encoding="utf-8")
    print("[SAVED]", OUT)
    print("[SAVED]", DISJOINT_OUT)
    print("[SAVED]", REPORT)
    print("[SAVED]", PLOT)
    if len(failed_disjoint):
        raise SystemExit(f"Cluster split disjointness check failed: {len(failed_disjoint)} rows. See {DISJOINT_OUT}")


if __name__ == "__main__":
    main()
