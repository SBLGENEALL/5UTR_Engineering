from __future__ import annotations

from pathlib import Path
from collections import Counter
from itertools import product
import argparse
import math
import shutil
import subprocess
import tempfile
import warnings

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import joblib

from scipy.stats import spearmanr, pearsonr
from sklearn.model_selection import train_test_split
from sklearn.pipeline import Pipeline
from sklearn.impute import SimpleImputer
from sklearn.ensemble import (
    ExtraTreesRegressor,
    RandomForestRegressor,
    HistGradientBoostingRegressor,
    ExtraTreesClassifier,
    RandomForestClassifier,
    HistGradientBoostingClassifier,
)
from sklearn.metrics import (
    r2_score,
    mean_squared_error,
    mean_absolute_error,
    roc_auc_score,
    average_precision_score,
)

warnings.filterwarnings("ignore")

BASE = Path.cwd()
SEQ = "utr5_sequence_tss_corrected"

LABEL_CANDIDATES = [
    BASE / "04_te_labeling/tables/tss_corrected_5utr_multiomics_labels.csv",
    BASE / "04_te_labeling/tables/tss_corrected_5utr_robust_public_te_labels.csv",
]

OUT_TABLE = BASE / "06_modeling/tables/heavy_rnafold_kmer6_model_search_results.csv"
OUT_WEIGHTS = BASE / "06_modeling/tables/heavy_rnafold_kmer6_selected_model_weights.csv"
OUT_CAND = BASE / "06_modeling/tables/heavy_rnafold_kmer6_50_100_candidate_scores.csv"
OUT_SUMMARY = BASE / "06_modeling/tables/heavy_rnafold_kmer6_summary.txt"
OUT_PLOT = BASE / "06_modeling/plots/heavy_rnafold_kmer6_gene_split_performance.png"
OUT_LIB = BASE / "07_library_design/tables/selected_1000_50_100bp_HEAVY_RNAfold_kmer6_library.csv"
OUT_FASTA = BASE / "07_library_design/fasta/selected_1000_50_100bp_HEAVY_RNAfold_kmer6_library.fasta"
OUT_LIB_SUMMARY = BASE / "07_library_design/qc/selected_1000_50_100bp_HEAVY_RNAfold_kmer6_summary.txt"
CACHE = BASE / "05_feature_extraction/rnafold/rnafold_features_cache.csv"

for p in [OUT_TABLE.parent, OUT_PLOT.parent, OUT_LIB.parent, OUT_FASTA.parent, OUT_LIB_SUMMARY.parent, CACHE.parent]:
    p.mkdir(parents=True, exist_ok=True)


def clean_seq(x):
    if pd.isna(x):
        return ""
    return str(x).upper().replace("U", "T").replace(" ", "").replace("\n", "")


def gc_content(seq):
    seq = clean_seq(seq)
    return (seq.count("G") + seq.count("C")) / len(seq) if len(seq) else np.nan


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
            for j in range(i + 3, len(seq) - 2, 3):
                if seq[j:j+3] in stops:
                    n += 1
                    longest = max(longest, j + 3 - i)
                    break
    return n, longest


def kmer_features(seq, ks=(3, 4, 5, 6), prefix=""):
    seq = clean_seq(seq)
    out = {}
    alphabet = "ACGT"
    for k in ks:
        all_kmers = ["".join(x) for x in product(alphabet, repeat=k)]
        cnt = Counter(seq[i:i+k] for i in range(max(0, len(seq)-k+1)) if set(seq[i:i+k]) <= set(alphabet))
        total = sum(cnt.values())
        for km in all_kmers:
            out[f"{prefix}{k}mer_{km}"] = cnt.get(km, 0) / total if total else 0.0
    return out


def onehot100(seq, prefix="onehot_"):
    seq = clean_seq(seq)[:100]
    out = {}
    for i in range(100):
        b = seq[i] if i < len(seq) else "N"
        for x in "ACGT":
            out[f"{prefix}{i:03d}_{x}"] = 1.0 if b == x else 0.0
    return out


def fold_single_chunk(records, exe="RNAfold"):
    """Run RNAfold on a chunk of [(id, seq)] and return {id: mfe}."""
    fasta_lines = []
    for rid, seq in records:
        fasta_lines.append(f">{rid}")
        fasta_lines.append(clean_seq(seq).replace("T", "U"))
    payload = "\n".join(fasta_lines) + "\n"

    r = subprocess.run(
        [exe, "--noPS"],
        input=payload,
        text=True,
        capture_output=True,
        check=True,
    )

    lines = [x.strip() for x in r.stdout.splitlines() if x.strip()]
    out = {}
    i = 0
    while i < len(lines):
        if not lines[i].startswith(">"):
            i += 1
            continue
        rid = lines[i][1:].split()[0]
        if i + 2 < len(lines):
            struct_line = lines[i + 2]
            # Typical: "....((((...)))) (-12.30)"
            mfe = np.nan
            if "(" in struct_line and ")" in struct_line:
                try:
                    mfe = float(struct_line.rsplit("(", 1)[1].split(")", 1)[0].strip())
                except Exception:
                    mfe = np.nan
            out[rid] = mfe
        i += 3
    return out


def rnafold_features_for_sequences(seqs, exe="RNAfold", chunk_size=500, max_rows=None):
    """Return dataframe: sequence, rnafold_mfe, rnafold_mfe_per_nt.

    Uses cache so reruns are fast.
    """
    seqs = [clean_seq(s) for s in seqs if clean_seq(s)]
    unique = list(dict.fromkeys(seqs))
    if max_rows is not None:
        unique = unique[:max_rows]

    if CACHE.exists():
        cache = pd.read_csv(CACHE)
    else:
        cache = pd.DataFrame(columns=["sequence", "rnafold_mfe", "rnafold_mfe_per_nt"])

    have = set(cache["sequence"].astype(str)) if len(cache) else set()
    todo = [s for s in unique if s not in have]

    if todo:
        if shutil.which(exe) is None:
            raise SystemExit(f"RNAfold executable not found: {exe}. Check PATH or conda env.")
        print(f"[RNAfold] cached={len(have):,}, new={len(todo):,}")
        rows = []
        for start in range(0, len(todo), chunk_size):
            chunk = todo[start:start + chunk_size]
            records = [(f"seq_{start+i}", s) for i, s in enumerate(chunk)]
            mfe_map = fold_single_chunk(records, exe=exe)
            for rid, seq in records:
                mfe = mfe_map.get(rid, np.nan)
                rows.append({
                    "sequence": seq,
                    "rnafold_mfe": mfe,
                    "rnafold_mfe_per_nt": mfe / len(seq) if len(seq) and pd.notna(mfe) else np.nan,
                })
            print(f"  RNAfold {min(start+chunk_size, len(todo)):,}/{len(todo):,}")
        new = pd.DataFrame(rows)
        cache = pd.concat([cache, new], ignore_index=True).drop_duplicates(subset=["sequence"], keep="last")
        cache.to_csv(CACHE, index=False)
        print("[SAVED RNAfold cache]", CACHE, cache.shape)

    return cache


def make_features(df, kmax=6, use_rnafold=True, use_onehot=True, rnafold_exe="RNAfold"):
    seqs = df[SEQ].map(clean_seq).tolist()
    fold_map = {}
    if use_rnafold:
        fold_df = rnafold_features_for_sequences(seqs, exe=rnafold_exe)
        fold_map = fold_df.set_index("sequence")[["rnafold_mfe", "rnafold_mfe_per_nt"]].to_dict("index")

    k_full = tuple(range(3, kmax + 1))
    k_tail = tuple(range(3, min(kmax, 5) + 1))

    rows = []
    for i, seq in enumerate(seqs):
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
            "aug_density": seq.count("ATG") / L if L else 0,
            "uorf_count": uorf,
            "longest_uorf_len": longest,
        }

        for b in "ACGT":
            r[f"{b.lower()}_fraction"] = seq.count(b) / L if L else 0.0
            r[f"tail100_{b.lower()}_fraction"] = tail100.count(b) / len(tail100) if len(tail100) else 0.0
            r[f"poly{b}_max_run"] = longest_run(seq, b)
            r[f"tail100_poly{b}_max_run"] = longest_run(tail100, b)

        for w in [20, 30, 50, 100]:
            head = seq[:w]
            tail = seq[-w:] if L >= w else seq
            r[f"head{w}_gc"] = gc_content(head)
            r[f"tail{w}_gc"] = gc_content(tail)
            r[f"head{w}_uaug"] = head.count("ATG")
            r[f"tail{w}_uaug"] = tail.count("ATG")

        r.update(kmer_features(seq, ks=k_full, prefix="full_"))
        r.update(kmer_features(tail100, ks=k_tail, prefix="tail100_"))
        r.update(kmer_features(head100, ks=(3, 4, 5), prefix="head100_"))

        if use_onehot:
            r.update(onehot100(tail100, prefix="tail100_onehot_"))

        if use_rnafold:
            fold = fold_map.get(seq, {})
            r["rnafold_mfe"] = fold.get("rnafold_mfe", np.nan)
            r["rnafold_mfe_per_nt"] = fold.get("rnafold_mfe_per_nt", np.nan)

        rows.append(r)
        if (i + 1) % 3000 == 0:
            print(f"  features {i+1:,}/{len(df):,}")

    return pd.DataFrame(rows)


def make_regressor(name, n_estimators=2000):
    if name == "ExtraTrees":
        return Pipeline([
            ("imputer", SimpleImputer(strategy="median")),
            ("model", ExtraTreesRegressor(
                n_estimators=n_estimators, max_features="sqrt",
                min_samples_leaf=3, random_state=42, n_jobs=-1
            )),
        ])
    if name == "RandomForest":
        return Pipeline([
            ("imputer", SimpleImputer(strategy="median")),
            ("model", RandomForestRegressor(
                n_estimators=n_estimators, max_features="sqrt",
                min_samples_leaf=3, random_state=42, n_jobs=-1
            )),
        ])
    if name == "HistGradientBoosting":
        return Pipeline([
            ("imputer", SimpleImputer(strategy="median")),
            ("model", HistGradientBoostingRegressor(
                max_iter=800, learning_rate=0.035, max_leaf_nodes=31,
                l2_regularization=0.02, random_state=42
            )),
        ])
    raise ValueError(name)


def make_classifier(name, n_estimators=2000):
    if name == "ExtraTrees":
        return Pipeline([
            ("imputer", SimpleImputer(strategy="median")),
            ("model", ExtraTreesClassifier(
                n_estimators=n_estimators, max_features="sqrt",
                min_samples_leaf=3, class_weight="balanced",
                random_state=42, n_jobs=-1
            )),
        ])
    if name == "RandomForest":
        return Pipeline([
            ("imputer", SimpleImputer(strategy="median")),
            ("model", RandomForestClassifier(
                n_estimators=n_estimators, max_features="sqrt",
                min_samples_leaf=3, class_weight="balanced",
                random_state=42, n_jobs=-1
            )),
        ])
    if name == "HistGradientBoosting":
        return Pipeline([
            ("imputer", SimpleImputer(strategy="median")),
            ("model", HistGradientBoostingClassifier(
                max_iter=800, learning_rate=0.035, max_leaf_nodes=31,
                l2_regularization=0.02, random_state=42
            )),
        ])
    raise ValueError(name)


def split_idx(df, mode):
    if mode == "random":
        idx = np.arange(len(df))
        return train_test_split(idx, test_size=0.2, random_state=42)
    groups = df["gene_name"].astype(str).fillna("NA").values if "gene_name" in df.columns else df["utr_id"].astype(str).values
    uniq = np.array(sorted(pd.unique(groups)))
    tr_g, te_g = train_test_split(uniq, test_size=0.2, random_state=42)
    return np.where(np.isin(groups, tr_g))[0], np.where(np.isin(groups, te_g))[0]


def reg_metrics(y, p):
    out = {
        "spearman": spearmanr(y, p).correlation if len(y) >= 3 else np.nan,
        "pearson": pearsonr(y, p)[0] if len(y) >= 3 else np.nan,
        "r2": r2_score(y, p),
        "rmse": math.sqrt(mean_squared_error(y, p)),
        "mae": mean_absolute_error(y, p),
    }
    for frac in [0.05, 0.10, 0.20]:
        k = max(1, int(len(y) * frac))
        top = np.argsort(p)[-k:]
        thr = np.quantile(y, 1 - frac)
        out[f"top{int(frac*100)}_enrichment"] = float(np.mean(y[top] >= thr))
    return out


def load_label():
    for p in LABEL_CANDIDATES:
        if p.exists():
            print("[LOAD]", p)
            return pd.read_csv(p), p
    raise SystemExit("No label table found. Run publicTE/multiomics pipeline first.")


def prep(df):
    df = df.copy()
    df[SEQ] = df[SEQ].map(clean_seq)
    if "utr5_length_final" in df.columns:
        df["length_for_filter"] = df["utr5_length_final"]
    elif "utr5_length_tss_corrected" in df.columns:
        df["length_for_filter"] = df["utr5_length_tss_corrected"]
    else:
        df["length_for_filter"] = df[SEQ].str.len()
    if "gc_content" not in df.columns:
        df["gc_content"] = df[SEQ].map(gc_content)
    if "uaug_count" not in df.columns:
        df["uaug_count"] = df[SEQ].str.count("ATG")
    if "has_proteomics_label" not in df.columns:
        df["has_proteomics_label"] = df.get("protein_abundance_rank", pd.Series(np.nan, index=df.index)).notna()
    return df


def filter_df(df, name):
    if name == "40_200_strict":
        mask = df["length_for_filter"].between(40, 200) & df["gc_content"].between(0.30, 0.75) & (df["uaug_count"] <= 1)
    elif name == "20_500_relaxed":
        mask = df["length_for_filter"].between(20, 500) & df["gc_content"].between(0.25, 0.80) & (df["uaug_count"] <= 3)
    elif name == "proteomics_20_500_relaxed":
        mask = (
            df["length_for_filter"].between(20, 500)
            & df["gc_content"].between(0.25, 0.80)
            & (df["uaug_count"] <= 3)
            & df["has_proteomics_label"].astype(bool)
        )
    elif name == "50_100_strict":
        mask = df["length_for_filter"].between(50, 100) & df["gc_content"].between(0.30, 0.75) & (df["uaug_count"] <= 1)
    else:
        raise ValueError(name)
    mask &= df[SEQ].str.len() > 0
    mask &= ~df[SEQ].str.contains("N", regex=False)
    return df[mask].copy().reset_index(drop=True)


def forbidden_sites(seq):
    sites = {
        "BsaI_GGTCTC": "GGTCTC",
        "BsaI_GAGACC": "GAGACC",
        "BsmBI_CGTCTC": "CGTCTC",
        "BsmBI_GAGACG": "GAGACG",
        "EcoRI_GAATTC": "GAATTC",
        "XhoI_CTCGAG": "CTCGAG",
        "NheI_GCTAGC": "GCTAGC",
        "AgeI_ACCGGT": "ACCGGT",
        "NotI_GCGGCCGC": "GCGGCCGC",
    }
    seq = clean_seq(seq)
    return ";".join([name for name, site in sites.items() if site in seq])


def write_fasta(df, path):
    with open(path, "w", encoding="utf-8") as out:
        for i, r in df.iterrows():
            uid = str(r.get("utr_id", f"row_{i}")).replace(" ", "_").replace("/", "_")
            gene = str(r.get("gene_name", "NA")).replace(" ", "_").replace("/", "_")
            group = str(r.get("heavy_group", "NA")).replace(" ", "_")
            seq = clean_seq(r[SEQ])
            out.write(f">{uid}|{gene}|len={len(seq)}|group={group}\n")
            for j in range(0, len(seq), 80):
                out.write(seq[j:j+80] + "\n")


def select_candidate_space(df):
    x = filter_df(df, "50_100_strict")
    x["forbidden_sites"] = x[SEQ].apply(forbidden_sites)
    x = x[
        x["forbidden_sites"].fillna("").astype(str).eq("")
        & x["robust_public_te_rank"].notna()
    ].copy()
    return x.drop_duplicates(subset=[SEQ]).reset_index(drop=True)


def rank_pct(s):
    return pd.to_numeric(s, errors="coerce").rank(pct=True)


def make_library(cand):
    df = cand.copy()
    for c in ["protein_abundance_rank", "protein_residual_rank", "multi_omics_utr_rank", "heavy_ensemble_score"]:
        if c not in df.columns:
            df[c] = np.nan

    df["protein_abundance_rank_fill"] = df["protein_abundance_rank"].fillna(df["robust_public_te_rank"])
    df["protein_residual_rank_fill"] = df["protein_residual_rank"].fillna(df["robust_public_te_rank"])
    df["multi_omics_utr_rank_fill"] = df["multi_omics_utr_rank"].fillna(df["robust_public_te_rank"])
    df["heavy_ensemble_score"] = df["heavy_ensemble_score"].fillna(df["robust_public_te_rank"])

    df["heavy_evidence_score"] = (
        0.42 * df["robust_public_te_rank"].fillna(0)
        + 0.18 * df["heavy_ensemble_score"].fillna(0)
        + 0.15 * df["protein_abundance_rank_fill"].fillna(0)
        + 0.12 * df["protein_residual_rank_fill"].fillna(0)
        + 0.08 * df.get("day_consensus_TE_rank", pd.Series(0, index=df.index)).fillna(0)
        + 0.05 * df.get("tss_confidence_score", pd.Series(0, index=df.index)).fillna(0)
    )

    selected = []
    used = set()

    def take(pool, n, group, sort_cols):
        nonlocal used
        x = pool[~pool[SEQ].isin(used)].copy()
        x = x.sort_values(sort_cols, ascending=[False] * len(sort_cols)).head(n)
        x["heavy_group"] = group
        selected.append(x)
        used.update(x[SEQ].tolist())

    take(df, 450, "A_heavy_evidence_top", ["heavy_evidence_score", "robust_public_te_rank"])
    take(df[df["robust_public_te_rank"].fillna(0) >= 0.80], 200, "B_publicTE_high_model_supported", ["heavy_ensemble_score", "robust_public_te_rank"])
    prot = df[df.get("has_proteomics_label", pd.Series(False, index=df.index)).astype(bool)].copy()
    take(prot, 150, "C_proteomics_supported", ["protein_abundance_rank_fill", "protein_residual_rank_fill", "heavy_evidence_score"])
    mid = df[df["robust_public_te_rank"].between(0.40, 0.85)].copy()
    take(mid, 175, "D_diverse_mid_high", ["heavy_ensemble_score", "heavy_evidence_score"])
    low = df[df["robust_public_te_rank"].fillna(1) <= 0.25].copy()
    take(low, 25, "E_low_TE_negative", ["gc_content"])

    out = pd.concat(selected).drop_duplicates(subset=[SEQ]) if selected else df.head(0)

    if len(out) < 1000:
        fill = df[~df[SEQ].isin(set(out[SEQ]))].sort_values(["heavy_evidence_score", "robust_public_te_rank"], ascending=False).head(1000 - len(out)).copy()
        fill["heavy_group"] = "F_fill_heavy_evidence"
        out = pd.concat([out, fill]).drop_duplicates(subset=[SEQ])

    out = out.head(1000).copy()
    out["library_index"] = range(1, len(out) + 1)
    front = [c for c in [
        "library_index", "heavy_group", "utr_id", "gene_id", "gene_name", SEQ,
        "length_for_filter", "gc_content", "uaug_count", "robust_public_te_rank",
        "multi_omics_utr_rank", "protein_abundance_rank", "protein_residual_rank",
        "heavy_ensemble_score", "heavy_evidence_score", "has_proteomics_label", "forbidden_sites"
    ] if c in out.columns]
    return out[front + [c for c in out.columns if c not in front]]


def plot_results(res):
    sub = res[res["split_mode"] == "gene_split"].copy()
    if sub.empty:
        return
    sub["label"] = sub["task"] + "\n" + sub["target"] + "\n" + sub["filter_name"] + "\n" + sub["model"]
    sub = sub.sort_values("selection_metric", ascending=True).tail(30)
    plt.figure(figsize=(14, max(8, 0.45 * len(sub))))
    plt.barh(sub["label"], sub["selection_metric"])
    plt.xlabel("Heavy selection metric")
    plt.title("Heavy RNAfold + k-mer6 model search")
    plt.tight_layout()
    plt.savefig(OUT_PLOT, dpi=220)
    plt.close()


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--n-estimators", type=int, default=2000)
    ap.add_argument("--kmax", type=int, default=6)
    ap.add_argument("--rnafold-exe", default="RNAfold")
    ap.add_argument("--no-rnafold", action="store_true")
    ap.add_argument("--no-onehot", action="store_true")
    ap.add_argument("--top-models", type=int, default=12)
    ap.add_argument("--filters", default="40_200_strict,20_500_relaxed,proteomics_20_500_relaxed")
    args = ap.parse_args()

    df, label_path = load_label()
    df = prep(df)

    filters = [x.strip() for x in args.filters.split(",") if x.strip()]
    targets = [c for c in ["robust_public_te_rank", "multi_omics_utr_rank", "protein_abundance_rank", "protein_residual_rank"] if c in df.columns]
    models = ["ExtraTrees", "RandomForest", "HistGradientBoosting"]
    rows = []
    feature_cache = {}

    for filt in filters:
        sub0 = filter_df(df, filt)
        print("\n" + "=" * 100)
        print("[FILTER]", filt, "rows", len(sub0), "proteomics", int(sub0.get("has_proteomics_label", pd.Series(False, index=sub0.index)).sum()))
        print("=" * 100)
        if len(sub0) < 500:
            print("[SKIP] too few rows")
            continue

        X0 = make_features(
            sub0,
            kmax=args.kmax,
            use_rnafold=not args.no_rnafold,
            use_onehot=not args.no_onehot,
            rnafold_exe=args.rnafold_exe,
        )
        feature_cache[filt] = (sub0, X0)

        for target in targets:
            valid = pd.to_numeric(sub0[target], errors="coerce").notna()
            sub = sub0[valid].copy().reset_index(drop=True)
            X = X0.loc[valid].reset_index(drop=True)
            if len(sub) < 500:
                continue
            if target.startswith("protein") and int(sub.get("has_proteomics_label", pd.Series(False, index=sub.index)).sum()) < 500:
                continue

            y = pd.to_numeric(sub[target], errors="coerce").values

            for split_mode in ["random", "gene_split"]:
                tr, te = split_idx(sub, split_mode)
                for model_name in models:
                    print("[REG]", filt, target, split_mode, model_name)
                    model = make_regressor(model_name, n_estimators=args.n_estimators)
                    model.fit(X.iloc[tr], y[tr])
                    pred = np.clip(model.predict(X.iloc[te]), 0, 1)
                    met = reg_metrics(y[te], pred)
                    metric = 0.50 * max(0, met["spearman"]) + 0.25 * met["top10_enrichment"] + 0.25 * met["top20_enrichment"]
                    rows.append({
                        "task": "regression", "filter_name": filt, "target": target,
                        "split_mode": split_mode, "model": model_name,
                        "n_estimators": args.n_estimators, "kmax": args.kmax,
                        "use_rnafold": not args.no_rnafold, "use_onehot": not args.no_onehot,
                        "n_rows": len(sub), "n_train": len(tr), "n_test": len(te),
                        "selection_metric": metric, **met
                    })
                    print("  metric", round(metric, 4), "spearman", round(met["spearman"], 4), "top10", round(met["top10_enrichment"], 4))

            # high/low classification
            hi = np.quantile(y, 0.80)
            lo = np.quantile(y, 0.30)
            cmask = (y >= hi) | (y <= lo)
            if cmask.sum() >= 300:
                cdf = sub.loc[cmask].copy().reset_index(drop=True)
                Xc = X.loc[cmask].reset_index(drop=True)
                yc = (pd.to_numeric(cdf[target], errors="coerce").values >= hi).astype(int)
                if len(np.unique(yc)) == 2:
                    for split_mode in ["random", "gene_split"]:
                        if split_mode == "random":
                            idx = np.arange(len(cdf))
                            tr, te = train_test_split(idx, test_size=0.2, random_state=42, stratify=yc)
                        else:
                            tr, te = split_idx(cdf, split_mode)
                        for model_name in models:
                            print("[CLF]", filt, target, split_mode, model_name)
                            clf = make_classifier(model_name, n_estimators=args.n_estimators)
                            clf.fit(Xc.iloc[tr], yc[tr])
                            prob = clf.predict_proba(Xc.iloc[te])[:, 1]
                            roc = roc_auc_score(yc[te], prob)
                            ap_score = average_precision_score(yc[te], prob)
                            metric = 0.45 * roc + 0.55 * ap_score
                            rows.append({
                                "task": "classification", "filter_name": filt, "target": target,
                                "split_mode": split_mode, "model": model_name,
                                "n_estimators": args.n_estimators, "kmax": args.kmax,
                                "use_rnafold": not args.no_rnafold, "use_onehot": not args.no_onehot,
                                "n_rows": len(cdf), "n_train": len(tr), "n_test": len(te),
                                "selection_metric": metric,
                                "roc_auc": roc, "average_precision": ap_score,
                                "spearman": np.nan, "pearson": np.nan, "r2": np.nan,
                                "rmse": np.nan, "mae": np.nan,
                                "top5_enrichment": np.nan, "top10_enrichment": np.nan, "top20_enrichment": np.nan,
                            })
                            print("  metric", round(metric, 4), "roc", round(roc, 4), "ap", round(ap_score, 4))

    res = pd.DataFrame(rows)
    res.to_csv(OUT_TABLE, index=False)
    plot_results(res)

    gene = res[res["split_mode"] == "gene_split"].sort_values("selection_metric", ascending=False).copy()
    selected = gene.head(args.top_models).copy()
    selected.to_csv(OUT_WEIGHTS, index=False)

    # Final candidate scoring
    cand = select_candidate_space(df)
    print("\n[FINAL CANDIDATES]", len(cand))
    pred_cols = []
    weights = []
    for i, row in selected.reset_index(drop=True).iterrows():
        filt = row["filter_name"]
        target = row["target"]
        task = row["task"]
        model_name = row["model"]

        train0 = filter_df(df, filt)
        valid = pd.to_numeric(train0[target], errors="coerce").notna()
        train = train0[valid].copy().reset_index(drop=True)
        if task == "classification":
            yraw = pd.to_numeric(train[target], errors="coerce").values
            hi = np.quantile(yraw, 0.80)
            lo = np.quantile(yraw, 0.30)
            cmask = (yraw >= hi) | (yraw <= lo)
            train = train.loc[cmask].copy().reset_index(drop=True)
            y = (pd.to_numeric(train[target], errors="coerce").values >= hi).astype(int)
            if len(np.unique(y)) < 2:
                continue
        else:
            y = pd.to_numeric(train[target], errors="coerce").values

        Xtrain = make_features(train, kmax=args.kmax, use_rnafold=not args.no_rnafold, use_onehot=not args.no_onehot, rnafold_exe=args.rnafold_exe)
        Xcand = make_features(cand, kmax=args.kmax, use_rnafold=not args.no_rnafold, use_onehot=not args.no_onehot, rnafold_exe=args.rnafold_exe)

        col = f"heavy_model_{i+1:02d}_{task}_{target}_{filt}_{model_name}"
        print("[FINAL FIT]", col)
        if task == "classification":
            m = make_classifier(model_name, n_estimators=args.n_estimators)
            m.fit(Xtrain, y)
            pred = m.predict_proba(Xcand)[:, 1]
        else:
            m = make_regressor(model_name, n_estimators=args.n_estimators)
            m.fit(Xtrain, y)
            pred = np.clip(m.predict(Xcand), 0, 1)

        cand[col] = pred
        pred_cols.append(col)
        weights.append(max(float(row["selection_metric"]), 0.0))
        joblib.dump({"model": m, "row": row.to_dict()}, BASE / f"06_modeling/models/{col}.joblib")

    if pred_cols:
        w = np.array(weights)
        if w.sum() <= 0:
            w = np.ones(len(pred_cols))
        w = w / w.sum()
        cand["heavy_ensemble_score"] = cand[pred_cols].fillna(0).values.dot(w)
    else:
        cand["heavy_ensemble_score"] = cand["robust_public_te_rank"]

    cand["heavy_ensemble_rank"] = rank_pct(cand["heavy_ensemble_score"])
    cand.to_csv(OUT_CAND, index=False)

    lib = make_library(cand)
    lib.to_csv(OUT_LIB, index=False)
    write_fasta(lib, OUT_FASTA)

    OUT_LIB_SUMMARY.write_text(
        "Heavy RNAfold kmer6 selected library summary\n"
        + "="*90 + "\n"
        + f"Total selected: {len(lib)}\n"
        + f"Candidate space: {len(cand)}\n\n"
        + "[Group counts]\n" + lib["heavy_group"].value_counts().to_string()
        + "\n\n[Robust public TE rank]\n" + lib["robust_public_te_rank"].describe().to_string()
        + "\n\n[Heavy ensemble score]\n" + lib["heavy_ensemble_score"].describe().to_string()
        + "\n\n[has proteomics label]\n" + (lib["has_proteomics_label"].value_counts(dropna=False).to_string() if "has_proteomics_label" in lib.columns else "NA")
        + "\n\n[uAUG]\n" + lib["uaug_count"].value_counts().sort_index().to_string()
        + "\n",
        encoding="utf-8"
    )

    OUT_SUMMARY.write_text(
        "Heavy RNAfold + kmer6 AutoML summary\n"
        + "="*90 + "\n"
        + f"label_path: {label_path}\n"
        + f"n_estimators: {args.n_estimators}\n"
        + f"kmax: {args.kmax}\n"
        + f"use_rnafold: {not args.no_rnafold}\n"
        + f"use_onehot: {not args.no_onehot}\n\n"
        + "[Selected gene-split models]\n"
        + selected.to_string(index=False)
        + "\n\n[Top 30 gene-split rows]\n"
        + gene.head(30).to_string(index=False)
        + f"\n\nSaved results: {OUT_TABLE}\nSaved weights: {OUT_WEIGHTS}\nSaved candidate scores: {OUT_CAND}\nSaved final library: {OUT_LIB}\nSaved plot: {OUT_PLOT}\n",
        encoding="utf-8"
    )

    print("[SAVED]", OUT_TABLE)
    print("[SAVED]", OUT_WEIGHTS)
    print("[SAVED]", OUT_CAND)
    print("[SAVED]", OUT_SUMMARY)
    print("[SAVED]", OUT_LIB)
    print("[SAVED]", OUT_FASTA)
    print("[SAVED]", OUT_LIB_SUMMARY)


if __name__ == "__main__":
    main()
