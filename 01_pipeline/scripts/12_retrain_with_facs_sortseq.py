from __future__ import annotations

"""Join the selected CHO 5'UTR library to Sort-seq phenotypes and retrain models.

The script deliberately keeps two experimental endpoints separate:

* high15_probability: primary high-expression-tail endpoint from SortSeq.
* expected_bin_score: secondary 6-to-1 whole-distribution endpoint.

Model performance is reported only from out-of-fold predictions.  Final models are
fit after cross-validation for later scoring, but their in-sample fit is never used
as the reported performance estimate.
"""

import argparse
from collections import Counter
from dataclasses import dataclass
from hashlib import sha256
from itertools import product
import json
import math
import os
from pathlib import Path
import re
import tempfile
from typing import Iterable

os.environ.setdefault(
    "MPLCONFIGDIR", str(Path(tempfile.gettempdir()) / "cho5utr_matplotlib_cache")
)

import joblib
import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy.stats import fisher_exact, kendalltau, pearsonr, spearmanr
from sklearn.ensemble import ExtraTreesRegressor
from sklearn.impute import SimpleImputer
from sklearn.metrics import (
    average_precision_score,
    mean_absolute_error,
    mean_squared_error,
    r2_score,
    roc_auc_score,
)
from sklearn.model_selection import GroupKFold
from sklearn.pipeline import Pipeline


SCRIPT_VERSION = "2.0.0"
BASE = Path.cwd()
DEFAULT_LIBRARY = (
    BASE
    / "07_library_design/tables/selected_2000_50_100bp_cluster_diverse_evidence_balanced_library.csv"
)
DEFAULT_OUTDIR = BASE / "08_reports/facs_retraining"

PRIMARY_TARGET = "high15_probability"
SECONDARY_TARGET = "expected_bin_score"
DEFAULT_TARGETS = [PRIMARY_TARGET, SECONDARY_TARGET]

SEQUENCE_COLUMN_CANDIDATES = [
    "utr5_sequence_tss_corrected",
    "sequence",
    "target_sequence",
    "sequence_key",
    "variant_id",
]

OLD_SCORE_CANDIDATES = [
    "heavy_ensemble_score",
    "model_support_score",
    "cluster_diverse_evidence_score",
    "robust_public_te_rank",
    "day_consensus_TE_rank",
    "multi_omics_utr_rank",
    "protein_abundance_rank",
    "protein_residual_rank",
    "automl_ensemble_score",
    "proteomics_enriched_score",
    "model_pred_rank_40_200train",
]

TRUE_VALUES = {"1", "true", "t", "yes", "y"}
FALSE_VALUES = {"0", "false", "f", "no", "n", ""}


def parse_args(argv: list[str] | None = None) -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description=(
            "Join the v1.4 selected 5'UTR library to SortSeq phenotypes, compare "
            "the original prediction scores with FACS outcomes, and train leakage-aware models."
        )
    )
    parser.add_argument("--selected-library", type=Path, default=DEFAULT_LIBRARY)
    parser.add_argument("--sortseq-results", type=Path, required=True)
    parser.add_argument("--outdir", type=Path, default=DEFAULT_OUTDIR)
    parser.add_argument("--library-sequence-column", default=None)
    parser.add_argument("--sortseq-sequence-column", default=None)
    parser.add_argument("--category-column", default="library_group")
    parser.add_argument("--eligibility-column", default="top15_read_support_pass")
    parser.add_argument("--min-total-bin-count", type=int, default=200)
    parser.add_argument("--min-match-rate", type=float, default=0.95)
    parser.add_argument("--expected-category-count", type=int, default=8)
    parser.add_argument(
        "--targets",
        default=",".join(DEFAULT_TARGETS),
        help="Comma-separated experimental targets; defaults to High15 and expected bin score.",
    )
    parser.add_argument(
        "--old-score-columns",
        default="",
        help="Extra comma-separated original score columns to compare and use in prior/hybrid models.",
    )
    parser.add_argument("--top-fractions", default="0.05,0.10,0.15,0.20")
    parser.add_argument("--cv-folds", type=int, default=5)
    parser.add_argument("--kmer-max", type=int, default=6)
    parser.add_argument("--n-estimators", type=int, default=300)
    parser.add_argument("--min-samples-leaf", type=int, default=3)
    parser.add_argument(
        "--n-jobs", type=int, default=max(1, min(32, os.cpu_count() or 1))
    )
    parser.add_argument("--random-seed", type=int, default=20260824)
    parser.add_argument("--skip-modeling", action="store_true")
    parser.add_argument("--no-final-models", action="store_true")
    parser.add_argument("--no-plots", action="store_true")
    return parser.parse_args(argv)


def read_table(path: Path) -> pd.DataFrame:
    if not path.exists():
        raise FileNotFoundError(path)
    suffixes = [item.lower() for item in path.suffixes]
    if ".tsv" in suffixes or path.suffix.lower() in {".txt"}:
        return pd.read_csv(path, sep="\t", low_memory=False)
    if ".csv" in suffixes:
        return pd.read_csv(path, low_memory=False)
    return pd.read_csv(path, sep=None, engine="python", low_memory=False)


def clean_sequence(value: object) -> str:
    if pd.isna(value):
        return ""
    return re.sub(r"[^ACGTN]", "", str(value).upper().replace("U", "T"))


def looks_like_sequence(series: pd.Series) -> bool:
    values = series.dropna().astype(str).head(200)
    if len(values) == 0:
        return False
    cleaned = values.map(clean_sequence)
    original_nonspace = values.str.replace(r"\s+", "", regex=True).str.len()
    return bool(((cleaned.str.len() >= 20) & (cleaned.str.len() == original_nonspace)).mean() >= 0.90)


def resolve_sequence_column(
    table: pd.DataFrame, requested: str | None, table_label: str
) -> str:
    if requested:
        if requested not in table.columns:
            raise ValueError(
                f"{table_label} sequence column {requested!r} is missing; "
                f"available columns include {list(table.columns[:20])}"
            )
        return requested
    for column in SEQUENCE_COLUMN_CANDIDATES:
        if column in table.columns and looks_like_sequence(table[column]):
            return column
    raise ValueError(
        f"Could not identify a sequence column in {table_label}. Supply an explicit column with "
        f"--{'library' if table_label == 'selected library' else 'sortseq'}-sequence-column."
    )


def to_bool(series: pd.Series, default: bool = False) -> pd.Series:
    if pd.api.types.is_bool_dtype(series):
        return series.fillna(default).astype(bool)
    normalized = series.fillna(str(default)).astype(str).str.strip().str.lower()
    unknown = ~normalized.isin(TRUE_VALUES | FALSE_VALUES)
    if bool(unknown.any()):
        bad = sorted(normalized[unknown].unique().tolist())[:5]
        raise ValueError(f"Unrecognized boolean values: {bad}")
    return normalized.isin(TRUE_VALUES)


@dataclass
class JoinResult:
    joined: pd.DataFrame
    unmatched_library: pd.DataFrame
    unmatched_sortseq: pd.DataFrame
    library_sequence_column: str
    sortseq_sequence_column: str
    match_rate: float


def join_library_to_sortseq(
    library: pd.DataFrame,
    sortseq: pd.DataFrame,
    library_sequence_column: str | None = None,
    sortseq_sequence_column: str | None = None,
    min_match_rate: float = 0.95,
) -> JoinResult:
    library = library.copy()
    sortseq = sortseq.copy()
    lib_seq = resolve_sequence_column(library, library_sequence_column, "selected library")
    facs_seq = resolve_sequence_column(sortseq, sortseq_sequence_column, "SortSeq results")

    library["_sequence_key"] = library[lib_seq].map(clean_sequence)
    sortseq["_sequence_key"] = sortseq[facs_seq].map(clean_sequence)
    if library["_sequence_key"].eq("").any():
        raise ValueError("The selected library contains empty/invalid sequence keys")
    if library["_sequence_key"].duplicated().any():
        examples = library.loc[
            library["_sequence_key"].duplicated(keep=False), "_sequence_key"
        ].head(5)
        raise ValueError(f"Selected-library sequence keys are not unique: {examples.tolist()}")
    sortseq_nonempty = sortseq[sortseq["_sequence_key"].ne("")].copy()
    if sortseq_nonempty["_sequence_key"].duplicated().any():
        examples = sortseq_nonempty.loc[
            sortseq_nonempty["_sequence_key"].duplicated(keep=False), "_sequence_key"
        ].head(5)
        raise ValueError(f"SortSeq sequence keys are not unique: {examples.tolist()}")

    merged = library.merge(
        sortseq_nonempty,
        on="_sequence_key",
        how="left",
        suffixes=("_selection", "_facs"),
        indicator=True,
        validate="one_to_one",
    )
    match_rate = float(merged["_merge"].eq("both").mean()) if len(merged) else 0.0
    unmatched_library = merged[merged["_merge"].ne("both")].copy()
    matched_keys = set(merged.loc[merged["_merge"].eq("both"), "_sequence_key"])
    unmatched_sortseq = sortseq_nonempty[
        ~sortseq_nonempty["_sequence_key"].isin(matched_keys)
    ].copy()
    if match_rate < min_match_rate:
        raise ValueError(
            f"Sequence match rate {match_rate:.3%} is below --min-match-rate "
            f"{min_match_rate:.3%}. Inspect unmatched_library.csv before modeling."
        )
    joined = merged[merged["_merge"].eq("both")].drop(columns="_merge").reset_index(drop=True)
    return JoinResult(
        joined=joined,
        unmatched_library=unmatched_library,
        unmatched_sortseq=unmatched_sortseq,
        library_sequence_column=lib_seq,
        sortseq_sequence_column=facs_seq,
        match_rate=match_rate,
    )


def resolve_merged_column(table: pd.DataFrame, requested: str) -> str:
    for candidate in [requested, f"{requested}_facs", f"{requested}_selection"]:
        if candidate in table.columns:
            return candidate
    raise ValueError(f"Required column {requested!r} is missing after the sequence join")


def resolve_selection_column(table: pd.DataFrame, requested: str) -> str:
    """Prefer the selected-library copy when both joined inputs contain a column."""

    for candidate in [requested, f"{requested}_selection", f"{requested}_facs"]:
        if candidate in table.columns:
            return candidate
    raise ValueError(f"Selected-library column {requested!r} is missing after the sequence join")


def add_analysis_annotations(
    joined: pd.DataFrame,
    targets: list[str],
    eligibility_column: str,
    min_total_bin_count: int,
) -> tuple[pd.DataFrame, pd.Series, dict[str, str]]:
    table = joined.copy()
    resolved_targets: dict[str, str] = {}
    for target in targets:
        column = resolve_merged_column(table, target)
        table[column] = pd.to_numeric(table[column], errors="coerce")
        resolved_targets[target] = column

    total_column = resolve_merged_column(table, "total_6bin_count")
    table[total_column] = pd.to_numeric(table[total_column], errors="coerce")
    eligible = table[total_column].ge(min_total_bin_count)
    eligibility_used = None
    for candidate in [
        eligibility_column,
        f"{eligibility_column}_facs",
        f"{eligibility_column}_selection",
    ]:
        if candidate in table.columns:
            eligible &= to_bool(table[candidate])
            eligibility_used = candidate
            break
    for column in resolved_targets.values():
        eligible &= table[column].notna()
    table["analysis_eligible"] = eligible
    table["analysis_eligibility_column"] = eligibility_used or "total_6bin_count_only"

    high15 = resolved_targets.get(PRIMARY_TARGET)
    expected = resolved_targets.get(SECONDARY_TARGET)
    if high15:
        table["facs_high15_percentile"] = np.nan
        table.loc[eligible, "facs_high15_percentile"] = (
            table.loc[eligible, high15].rank(method="average", pct=True) * 100
        )
        p = table["facs_high15_percentile"]
        table["facs_high15_tier"] = np.select(
            [p >= 95, p >= 85, p >= 75, p >= 50],
            ["Top 5%", "Top 5-15%", "Top 15-25%", "Top 25-50%"],
            default="Bottom 50%",
        )
        table.loc[~eligible, "facs_high15_tier"] = "Ineligible"
    if expected:
        table["facs_expected_score_percentile"] = np.nan
        table.loc[eligible, "facs_expected_score_percentile"] = (
            table.loc[eligible, expected].rank(method="average", pct=True) * 100
        )
    return table, eligible, resolved_targets


def choose_old_score_columns(
    table: pd.DataFrame, requested: Iterable[str] = (), min_non_null: int = 20
) -> list[str]:
    chosen: list[str] = []
    for base in [*OLD_SCORE_CANDIDATES, *requested]:
        for candidate in [base, f"{base}_selection"]:
            if candidate in table.columns:
                numeric = pd.to_numeric(table[candidate], errors="coerce")
                if int(numeric.notna().sum()) >= min_non_null and candidate not in chosen:
                    chosen.append(candidate)
                break
    return chosen


def benjamini_hochberg(pvalues: pd.Series) -> pd.Series:
    values = pd.to_numeric(pvalues, errors="coerce").to_numpy(dtype=float)
    result = np.full(len(values), np.nan)
    valid = np.flatnonzero(np.isfinite(values))
    if len(valid) == 0:
        return pd.Series(result, index=pvalues.index)
    order = valid[np.argsort(values[valid])]
    ranked = values[order] * len(order) / np.arange(1, len(order) + 1)
    ranked = np.minimum.accumulate(ranked[::-1])[::-1]
    result[order] = np.clip(ranked, 0, 1)
    return pd.Series(result, index=pvalues.index)


def safe_correlation(x: pd.Series, y: pd.Series, method: str) -> tuple[float, float]:
    pair = pd.concat(
        [pd.to_numeric(x, errors="coerce"), pd.to_numeric(y, errors="coerce")], axis=1
    ).dropna()
    if len(pair) < 3 or pair.iloc[:, 0].nunique() < 2 or pair.iloc[:, 1].nunique() < 2:
        return float("nan"), float("nan")
    if method == "spearman":
        result = spearmanr(pair.iloc[:, 0], pair.iloc[:, 1])
    elif method == "pearson":
        result = pearsonr(pair.iloc[:, 0], pair.iloc[:, 1])
    elif method == "kendall":
        result = kendalltau(pair.iloc[:, 0], pair.iloc[:, 1])
    else:
        raise ValueError(method)
    return float(result.statistic), float(result.pvalue)


def score_correlations(
    table: pd.DataFrame,
    eligible: pd.Series,
    old_scores: list[str],
    resolved_targets: dict[str, str],
    category_column: str,
) -> tuple[pd.DataFrame, pd.DataFrame]:
    rows: list[dict[str, object]] = []
    category_rows: list[dict[str, object]] = []
    for score in old_scores:
        old = pd.to_numeric(table[score], errors="coerce")
        for target_name, target_column in resolved_targets.items():
            mask = eligible & old.notna() & table[target_column].notna()
            x = old[mask]
            y = table.loc[mask, target_column]
            spearman, spearman_p = safe_correlation(x, y, "spearman")
            pearson, pearson_p = safe_correlation(x, y, "pearson")
            kendall, kendall_p = safe_correlation(x, y, "kendall")
            old_pct = x.rank(method="average", pct=True)
            target_pct = y.rank(method="average", pct=True)
            rows.append(
                {
                    "old_score": score,
                    "facs_target": target_name,
                    "n": int(mask.sum()),
                    "spearman_rho": spearman,
                    "spearman_p": spearman_p,
                    "pearson_r": pearson,
                    "pearson_p": pearson_p,
                    "kendall_tau": kendall,
                    "kendall_p": kendall_p,
                    "mean_absolute_percentile_gap": float(
                        (old_pct - target_pct).abs().mean() * 100
                    ),
                }
            )
            if category_column in table.columns:
                for category, indices in table[mask].groupby(category_column, dropna=False).groups.items():
                    if len(indices) < 10:
                        continue
                    rho, pvalue = safe_correlation(old.loc[indices], table.loc[indices, target_column], "spearman")
                    category_rows.append(
                        {
                            "category": category,
                            "old_score": score,
                            "facs_target": target_name,
                            "n": len(indices),
                            "spearman_rho": rho,
                            "spearman_p": pvalue,
                        }
                    )
    result = pd.DataFrame(rows)
    if len(result):
        result["spearman_fdr"] = benjamini_hochberg(result["spearman_p"])
        result["pearson_fdr"] = benjamini_hochberg(result["pearson_p"])
        result["kendall_fdr"] = benjamini_hochberg(result["kendall_p"])
    category_result = pd.DataFrame(category_rows)
    if len(category_result):
        category_result["spearman_fdr"] = benjamini_hochberg(
            category_result["spearman_p"]
        )
    return result, category_result


def topn_overlap(
    table: pd.DataFrame,
    eligible: pd.Series,
    old_scores: list[str],
    resolved_targets: dict[str, str],
    top_ns: Iterable[int] = (20, 50, 100, 200),
) -> pd.DataFrame:
    rows: list[dict[str, object]] = []
    for score in old_scores:
        old = pd.to_numeric(table[score], errors="coerce")
        for target_name, target_column in resolved_targets.items():
            mask = eligible & old.notna() & table[target_column].notna()
            valid = table.loc[mask].copy()
            valid["_old"] = old[mask]
            for requested_n in top_ns:
                n = min(int(requested_n), len(valid))
                if n < 1:
                    continue
                old_top = set(valid.nlargest(n, "_old").index)
                observed_top = set(valid.nlargest(n, target_column).index)
                overlap = len(old_top & observed_top)
                rows.append(
                    {
                        "old_score": score,
                        "facs_target": target_name,
                        "top_n": n,
                        "overlap_n": overlap,
                        "overlap_fraction": overlap / n,
                        "jaccard": overlap / len(old_top | observed_top),
                        "random_expected_overlap": n * n / len(valid),
                        "eligible_n": len(valid),
                    }
                )
    return pd.DataFrame(rows)


def category_tables(
    table: pd.DataFrame,
    eligible: pd.Series,
    category_column: str,
    resolved_targets: dict[str, str],
) -> tuple[pd.DataFrame, pd.DataFrame]:
    if category_column not in table.columns:
        raise ValueError(f"Category column {category_column!r} is missing")
    high15 = resolved_targets.get(PRIMARY_TARGET)
    expected = resolved_targets.get(SECONDARY_TARGET)
    rows: list[dict[str, object]] = []
    for category, group in table.groupby(category_column, dropna=False):
        group_eligible = group["analysis_eligible"]
        row: dict[str, object] = {
            "category": category,
            "selected_n": len(group),
            "selected_fraction": len(group) / len(table) if len(table) else np.nan,
            "eligible_n": int(group_eligible.sum()),
            "eligible_fraction_within_category": float(group_eligible.mean()),
        }
        eligible_group = group[group_eligible]
        if high15:
            row.update(
                {
                    "high15_mean": float(eligible_group[high15].mean()),
                    "high15_median": float(eligible_group[high15].median()),
                    "high15_q25": float(eligible_group[high15].quantile(0.25)),
                    "high15_q75": float(eligible_group[high15].quantile(0.75)),
                    "high15_top15_n": int(
                        eligible_group["facs_high15_percentile"].ge(85).sum()
                    ),
                    "high15_top15_rate": float(
                        eligible_group["facs_high15_percentile"].ge(85).mean()
                    )
                    if len(eligible_group)
                    else np.nan,
                }
            )
        if expected:
            row.update(
                {
                    "expected_score_mean": float(eligible_group[expected].mean()),
                    "expected_score_median": float(eligible_group[expected].median()),
                    "expected_score_q25": float(eligible_group[expected].quantile(0.25)),
                    "expected_score_q75": float(eligible_group[expected].quantile(0.75)),
                }
            )
        rows.append(row)
    summary = pd.DataFrame(rows).sort_values("selected_n", ascending=False)

    tiers = ["Top 5%", "Top 5-15%", "Top 15-25%", "Top 25-50%", "Bottom 50%"]
    tier_rows: list[dict[str, object]] = []
    if "facs_high15_tier" in table.columns:
        eligible_table = table[eligible]
        for category, group in eligible_table.groupby(category_column, dropna=False):
            counts = group["facs_high15_tier"].value_counts()
            for tier in tiers:
                count = int(counts.get(tier, 0))
                tier_rows.append(
                    {
                        "category": category,
                        "facs_high15_tier": tier,
                        "n": count,
                        "fraction_within_category": count / len(group) if len(group) else np.nan,
                    }
                )
    return summary, pd.DataFrame(tier_rows)


def category_hit_enrichment(
    table: pd.DataFrame,
    eligible: pd.Series,
    category_column: str,
    target_column: str,
    fractions: Iterable[float],
) -> pd.DataFrame:
    valid = table[eligible & table[target_column].notna()].copy()
    rows: list[dict[str, object]] = []
    if len(valid) == 0:
        return pd.DataFrame()
    categories = valid[category_column].fillna("NA")
    for fraction in fractions:
        hit_n = max(1, int(math.ceil(len(valid) * float(fraction))))
        hit_indices = set(valid.nlargest(hit_n, target_column).index)
        global_rate = hit_n / len(valid)
        for category, indices in categories.groupby(categories).groups.items():
            index_set = set(indices)
            group_n = len(index_set)
            observed = len(index_set & hit_indices)
            other_hits = hit_n - observed
            other_nonhits = (len(valid) - group_n) - other_hits
            odds_ratio, pvalue = fisher_exact(
                [[observed, group_n - observed], [other_hits, other_nonhits]],
                alternative="greater",
            )
            observed_rate = observed / group_n if group_n else np.nan
            rows.append(
                {
                    "top_fraction": float(fraction),
                    "category": category,
                    "eligible_n": group_n,
                    "observed_hits": observed,
                    "expected_hits": group_n * global_rate,
                    "hit_rate_within_category": observed_rate,
                    "global_hit_rate": global_rate,
                    "fold_enrichment": observed_rate / global_rate if global_rate else np.nan,
                    "odds_ratio": float(odds_ratio),
                    "fisher_greater_p": float(pvalue),
                }
            )
    result = pd.DataFrame(rows)
    if len(result):
        result["fisher_greater_fdr"] = benjamini_hochberg(result["fisher_greater_p"])
    return result


class DisjointSet:
    def __init__(self, n: int):
        self.parent = list(range(n))
        self.rank = [0] * n

    def find(self, item: int) -> int:
        while self.parent[item] != item:
            self.parent[item] = self.parent[self.parent[item]]
            item = self.parent[item]
        return item

    def union(self, left: int, right: int) -> None:
        left_root, right_root = self.find(left), self.find(right)
        if left_root == right_root:
            return
        if self.rank[left_root] < self.rank[right_root]:
            left_root, right_root = right_root, left_root
        self.parent[right_root] = left_root
        if self.rank[left_root] == self.rank[right_root]:
            self.rank[left_root] += 1


def choose_gene_column(table: pd.DataFrame) -> str | None:
    for base in ["gene_id", "gene_name"]:
        for candidate in [base, f"{base}_selection"]:
            if candidate in table.columns:
                values = table[candidate].dropna().astype(str).str.strip()
                if bool(values.ne("").any()):
                    return candidate
    return None


def choose_cluster_column(table: pd.DataFrame) -> str | None:
    for candidate in ["seq_cluster_id", "seq_cluster_id_selection"]:
        if candidate in table.columns:
            return candidate
    return None


def build_leakage_groups(table: pd.DataFrame) -> tuple[pd.Series, dict[str, object]]:
    """Create connected components so neither gene nor sequence cluster crosses folds."""

    reset = table.reset_index(drop=True)
    dsu = DisjointSet(len(reset))
    columns = [column for column in [choose_cluster_column(reset), choose_gene_column(reset)] if column]
    for column in columns:
        groups: dict[str, int] = {}
        for position, value in enumerate(reset[column]):
            if pd.isna(value):
                continue
            key = str(value).strip()
            if not key or key.lower() == "nan":
                continue
            if key in groups:
                dsu.union(position, groups[key])
            else:
                groups[key] = position
    roots = [dsu.find(position) for position in range(len(reset))]
    root_to_label = {root: f"component_{number:05d}" for number, root in enumerate(sorted(set(roots)), 1)}
    labels = pd.Series([root_to_label[root] for root in roots], index=table.index, name="cv_group")
    return labels, {
        "group_columns": columns,
        "n_components": int(labels.nunique()),
        "largest_component_n": int(labels.value_counts().max()) if len(labels) else 0,
    }


def sequence_feature_names(kmer_max: int) -> list[str]:
    if not 2 <= kmer_max <= 6:
        raise ValueError("--kmer-max must be between 2 and 6")
    names = [
        "length",
        "gc_content",
        "fraction_A",
        "fraction_C",
        "fraction_G",
        "fraction_T",
        "purine_fraction",
        "pyrimidine_fraction",
        "gc_first10",
        "gc_last10",
        "gc_first20",
        "gc_last20",
        "atg_per_100nt",
        "stop_per_100nt",
        "uorf_count",
        "longest_uorf_nt",
        "max_run_A",
        "max_run_C",
        "max_run_G",
        "max_run_T",
        "shannon_entropy",
    ]
    for k in range(2, kmer_max + 1):
        names.extend(f"k{k}_{''.join(chars)}" for chars in product("ACGT", repeat=k))
    return names


def longest_run(sequence: str, base: str) -> int:
    best = current = 0
    for letter in sequence:
        if letter == base:
            current += 1
            best = max(best, current)
        else:
            current = 0
    return best


def uorf_features(sequence: str) -> tuple[int, int]:
    stops = {"TAA", "TAG", "TGA"}
    count = longest = 0
    for start in range(max(0, len(sequence) - 2)):
        if sequence[start : start + 3] != "ATG":
            continue
        for stop in range(start + 3, len(sequence) - 2, 3):
            if sequence[stop : stop + 3] in stops:
                count += 1
                longest = max(longest, stop + 3 - start)
                break
    return count, longest


def gc_fraction(sequence: str) -> float:
    return (
        (sequence.count("G") + sequence.count("C")) / len(sequence)
        if sequence
        else 0.0
    )


def extract_sequence_features(sequences: pd.Series, kmer_max: int = 6) -> tuple[np.ndarray, list[str]]:
    names = sequence_feature_names(kmer_max)
    matrix = np.zeros((len(sequences), len(names)), dtype=np.float32)
    kmer_blocks: dict[int, tuple[int, dict[str, int]]] = {}
    offset = 21
    for k in range(2, kmer_max + 1):
        kmers = ["".join(chars) for chars in product("ACGT", repeat=k)]
        kmer_blocks[k] = (offset, {kmer: index for index, kmer in enumerate(kmers)})
        offset += len(kmers)

    for row_number, raw_sequence in enumerate(sequences):
        sequence = clean_sequence(raw_sequence).replace("N", "")
        length = len(sequence)
        if length == 0:
            continue
        counts = Counter(sequence)
        uorf_count, longest_uorf = uorf_features(sequence)
        probabilities = np.array([counts[base] / length for base in "ACGT"], dtype=float)
        nonzero = probabilities[probabilities > 0]
        entropy = float(-(nonzero * np.log2(nonzero)).sum()) if len(nonzero) else 0.0
        base_values = [
            length,
            gc_fraction(sequence),
            counts["A"] / length,
            counts["C"] / length,
            counts["G"] / length,
            counts["T"] / length,
            (counts["A"] + counts["G"]) / length,
            (counts["C"] + counts["T"]) / length,
            gc_fraction(sequence[:10]),
            gc_fraction(sequence[-10:]),
            gc_fraction(sequence[:20]),
            gc_fraction(sequence[-20:]),
            sequence.count("ATG") * 100 / length,
            sum(sequence.count(stop) for stop in ["TAA", "TAG", "TGA"]) * 100 / length,
            uorf_count,
            longest_uorf,
            longest_run(sequence, "A"),
            longest_run(sequence, "C"),
            longest_run(sequence, "G"),
            longest_run(sequence, "T"),
            entropy,
        ]
        matrix[row_number, :21] = np.asarray(base_values, dtype=np.float32)
        for k, (block_start, lookup) in kmer_blocks.items():
            denominator = max(0, length - k + 1)
            if denominator == 0:
                continue
            for position in range(denominator):
                kmer = sequence[position : position + k]
                index = lookup.get(kmer)
                if index is not None:
                    matrix[row_number, block_start + index] += 1.0 / denominator
    return matrix, names


def estimator_pipeline(
    seed: int, n_estimators: int, min_samples_leaf: int, n_jobs: int
) -> Pipeline:
    return Pipeline(
        [
            ("imputer", SimpleImputer(strategy="median")),
            (
                "estimator",
                ExtraTreesRegressor(
                    n_estimators=n_estimators,
                    min_samples_leaf=min_samples_leaf,
                    max_features="sqrt",
                    random_state=seed,
                    n_jobs=n_jobs,
                ),
            ),
        ]
    )


def regression_metrics(observed: pd.Series, predicted: np.ndarray) -> dict[str, float]:
    observed_array = pd.to_numeric(observed, errors="coerce").to_numpy(dtype=float)
    predicted_array = np.asarray(predicted, dtype=float)
    rho, _ = safe_correlation(pd.Series(observed_array), pd.Series(predicted_array), "spearman")
    pearson, _ = safe_correlation(pd.Series(observed_array), pd.Series(predicted_array), "pearson")
    top_binary = pd.Series(observed_array).rank(method="first", pct=True).ge(0.85).astype(int)
    metrics = {
        "spearman_rho": rho,
        "pearson_r": pearson,
        "mae": float(mean_absolute_error(observed_array, predicted_array)),
        "rmse": float(math.sqrt(mean_squared_error(observed_array, predicted_array))),
        "r2": float(r2_score(observed_array, predicted_array)),
    }
    if top_binary.nunique() == 2:
        metrics["top15_roc_auc"] = float(roc_auc_score(top_binary, predicted_array))
        metrics["top15_average_precision"] = float(
            average_precision_score(top_binary, predicted_array)
        )
    else:
        metrics["top15_roc_auc"] = np.nan
        metrics["top15_average_precision"] = np.nan
    for top_n in [20, 50, 100, 200]:
        n = min(top_n, len(observed_array))
        observed_top = set(np.argsort(observed_array)[-n:])
        predicted_top = set(np.argsort(predicted_array)[-n:])
        metrics[f"top{top_n}_overlap_fraction"] = len(observed_top & predicted_top) / n
    return metrics


@dataclass
class ModelingResult:
    oof_predictions: pd.DataFrame
    metrics: pd.DataFrame
    category_metrics: pd.DataFrame
    group_summary: dict[str, object]
    final_model_paths: list[Path]


def run_oof_models(
    table: pd.DataFrame,
    eligible: pd.Series,
    resolved_targets: dict[str, str],
    old_score_columns: list[str],
    sequence_column: str,
    category_column: str,
    outdir: Path,
    kmer_max: int = 6,
    cv_folds: int = 5,
    n_estimators: int = 300,
    min_samples_leaf: int = 3,
    n_jobs: int = 1,
    random_seed: int = 20260824,
    save_final_models: bool = True,
) -> ModelingResult:
    frame = table[eligible].copy().reset_index(drop=False).rename(columns={"index": "joined_row"})
    if len(frame) < max(30, cv_folds * 5):
        raise ValueError(f"Only {len(frame)} eligible rows; too few for {cv_folds}-fold modeling")
    groups, group_summary = build_leakage_groups(frame)
    n_groups = int(groups.nunique())
    n_splits = min(cv_folds, n_groups)
    if n_splits < 2:
        raise ValueError("Fewer than two independent gene/sequence-cluster groups")
    group_summary["requested_folds"] = cv_folds
    group_summary["used_folds"] = n_splits

    sequence_features, sequence_feature_columns = extract_sequence_features(
        frame[sequence_column], kmer_max=kmer_max
    )
    feature_sets: dict[str, tuple[np.ndarray, list[str]]] = {
        "sequence_only": (sequence_features, sequence_feature_columns)
    }
    usable_prior_columns = [
        column
        for column in old_score_columns
        if pd.to_numeric(frame[column], errors="coerce").notna().sum() >= 20
    ]
    if usable_prior_columns:
        prior = frame[usable_prior_columns].apply(pd.to_numeric, errors="coerce").to_numpy(dtype=float)
        feature_sets["prior_only"] = (prior, usable_prior_columns)
        feature_sets["hybrid"] = (
            np.hstack([sequence_features, prior]),
            [*sequence_feature_columns, *[f"prior::{column}" for column in usable_prior_columns]],
        )

    splitter = GroupKFold(n_splits=n_splits)
    split_indices = list(splitter.split(frame, groups=groups))
    fold_assignment = np.full(len(frame), -1, dtype=int)
    for fold, (_, test_indices) in enumerate(split_indices, start=1):
        fold_assignment[test_indices] = fold
    oof = frame[
        [
            column
            for column in [
                "joined_row",
                "library_index",
                "utr_id",
                category_column,
                sequence_column,
                "_sequence_key",
            ]
            if column in frame.columns
        ]
    ].copy()
    oof["cv_group"] = groups.to_numpy()
    oof["cv_fold"] = fold_assignment

    metric_rows: list[dict[str, object]] = []
    category_rows: list[dict[str, object]] = []
    final_model_paths: list[Path] = []
    model_dir = outdir / "models"
    if save_final_models:
        model_dir.mkdir(parents=True, exist_ok=True)

    for target_name, target_column in resolved_targets.items():
        y = pd.to_numeric(frame[target_column], errors="raise").reset_index(drop=True)
        oof[f"observed__{target_name}"] = y
        for model_name, (features, feature_names) in feature_sets.items():
            predictions = np.full(len(frame), np.nan, dtype=float)
            for fold, (train_indices, test_indices) in enumerate(split_indices, start=1):
                model = estimator_pipeline(
                    seed=random_seed + fold,
                    n_estimators=n_estimators,
                    min_samples_leaf=min_samples_leaf,
                    n_jobs=n_jobs,
                )
                model.fit(features[train_indices], y.iloc[train_indices])
                predictions[test_indices] = model.predict(features[test_indices])
            prediction_column = f"oof__{target_name}__{model_name}"
            oof[prediction_column] = predictions
            metrics = regression_metrics(y, predictions)
            metric_rows.append(
                {
                    "target": target_name,
                    "model": model_name,
                    "n": len(frame),
                    "n_features": features.shape[1],
                    "n_folds": n_splits,
                    **metrics,
                }
            )
            if category_column in frame.columns:
                for category, indices in frame.groupby(category_column, dropna=False).groups.items():
                    if len(indices) < 10:
                        continue
                    category_metric = regression_metrics(y.loc[indices], predictions[list(indices)])
                    category_rows.append(
                        {
                            "target": target_name,
                            "model": model_name,
                            "category": category,
                            "n": len(indices),
                            **category_metric,
                        }
                    )
            if save_final_models:
                final_model = estimator_pipeline(
                    seed=random_seed,
                    n_estimators=n_estimators,
                    min_samples_leaf=min_samples_leaf,
                    n_jobs=n_jobs,
                )
                final_model.fit(features, y)
                package = {
                    "script_version": SCRIPT_VERSION,
                    "target": target_name,
                    "model_name": model_name,
                    "model": final_model,
                    "feature_names": feature_names,
                    "kmer_max": kmer_max,
                    "old_score_columns": usable_prior_columns,
                    "sequence_column": sequence_column,
                    "training_n": len(frame),
                    "reported_performance_source": "group_disjoint_out_of_fold_predictions",
                }
                model_path = model_dir / f"facs_{target_name}__{model_name}.joblib"
                joblib.dump(package, model_path, compress=3)
                final_model_paths.append(model_path)

    metrics = pd.DataFrame(metric_rows).sort_values(
        ["target", "spearman_rho"], ascending=[True, False]
    )
    category_metrics = pd.DataFrame(category_rows)
    return ModelingResult(
        oof_predictions=oof,
        metrics=metrics,
        category_metrics=category_metrics,
        group_summary=group_summary,
        final_model_paths=final_model_paths,
    )


def write_plots(
    outdir: Path,
    table: pd.DataFrame,
    eligible: pd.Series,
    category_column: str,
    resolved_targets: dict[str, str],
    correlations: pd.DataFrame,
    hit_enrichment: pd.DataFrame,
    modeling: ModelingResult | None,
) -> list[Path]:
    plot_dir = outdir / "plots"
    plot_dir.mkdir(parents=True, exist_ok=True)
    written: list[Path] = []
    valid = table[eligible].copy()
    high15 = resolved_targets.get(PRIMARY_TARGET)
    expected = resolved_targets.get(SECONDARY_TARGET)

    counts = table[category_column].fillna("NA").value_counts().sort_values()
    fig, ax = plt.subplots(figsize=(9, max(4.5, 0.45 * len(counts))))
    ax.barh(counts.index.astype(str), counts.values, color="#4C78A8")
    for y, value in enumerate(counts.values):
        ax.text(value, y, f" {value}", va="center", fontsize=9)
    ax.set_xlabel("Selected UTR count")
    ax.set_title("v1.4 selected-library distribution by exact library_group")
    fig.tight_layout()
    path = plot_dir / "01_category_counts.png"
    fig.savefig(path, dpi=180)
    plt.close(fig)
    written.append(path)

    if high15:
        categories = counts.index.astype(str).tolist()
        data = [
            valid.loc[valid[category_column].fillna("NA").astype(str).eq(category), high15].dropna()
            for category in categories
        ]
        fig, ax = plt.subplots(figsize=(10, max(5, 0.55 * len(categories))))
        ax.boxplot(data, labels=categories, vert=False, showfliers=False)
        ax.set_xlabel("FACS High15 probability")
        ax.set_title("Experimental High15 distribution within each selection category")
        fig.tight_layout()
        path = plot_dir / "02_high15_by_category.png"
        fig.savefig(path, dpi=180)
        plt.close(fig)
        written.append(path)

    if high15 and expected:
        fig, ax = plt.subplots(figsize=(6.5, 5.4))
        ax.scatter(valid[expected], valid[high15], s=14, alpha=0.45, color="#2A9D8F")
        rho, _ = safe_correlation(valid[expected], valid[high15], "spearman")
        ax.set_xlabel("Expected bin score (6 to 1)")
        ax.set_ylabel("High15 probability")
        ax.set_title(f"Two experimental endpoints (Spearman rho={rho:.3f})")
        fig.tight_layout()
        path = plot_dir / "03_high15_vs_expected_score.png"
        fig.savefig(path, dpi=180)
        plt.close(fig)
        written.append(path)

    if high15 and len(correlations):
        ranked = correlations[correlations["facs_target"].eq(PRIMARY_TARGET)].copy()
        ranked = ranked.sort_values("spearman_rho", key=lambda x: x.abs(), ascending=False).head(6)
        if len(ranked):
            ncols = 2
            nrows = math.ceil(len(ranked) / ncols)
            fig, axes = plt.subplots(nrows, ncols, figsize=(10, 4.3 * nrows), squeeze=False)
            for ax, (_, row) in zip(axes.ravel(), ranked.iterrows()):
                score = row["old_score"]
                mask = eligible & pd.to_numeric(table[score], errors="coerce").notna()
                ax.scatter(table.loc[mask, score], table.loc[mask, high15], s=12, alpha=0.4)
                ax.set_xlabel(score)
                ax.set_ylabel("High15 probability")
                ax.set_title(f"rho={row['spearman_rho']:.3f}; n={int(row['n'])}")
            for ax in axes.ravel()[len(ranked) :]:
                ax.axis("off")
            fig.suptitle("Original scores versus experimental High15", y=1.01)
            fig.tight_layout()
            path = plot_dir / "04_original_scores_vs_high15.png"
            fig.savefig(path, dpi=180, bbox_inches="tight")
            plt.close(fig)
            written.append(path)

    if len(hit_enrichment):
        top15_rows = hit_enrichment[np.isclose(hit_enrichment["top_fraction"], 0.15)].copy()
        top15_rows = top15_rows.sort_values("fold_enrichment")
        fig, ax = plt.subplots(figsize=(9, max(4.5, 0.45 * len(top15_rows))))
        ax.barh(top15_rows["category"].astype(str), top15_rows["fold_enrichment"], color="#F28E2B")
        ax.axvline(1.0, color="black", linestyle="--", linewidth=1)
        ax.set_xlabel("Observed / expected Top15 hit rate")
        ax.set_title("Experimental Top15 enrichment by selection category")
        fig.tight_layout()
        path = plot_dir / "05_category_top15_enrichment.png"
        fig.savefig(path, dpi=180)
        plt.close(fig)
        written.append(path)

    if modeling is not None and len(modeling.metrics):
        best = modeling.metrics.sort_values("spearman_rho", ascending=False).iloc[0]
        target = str(best["target"])
        model = str(best["model"])
        observed_col = f"observed__{target}"
        predicted_col = f"oof__{target}__{model}"
        fig, ax = plt.subplots(figsize=(6.3, 5.4))
        ax.scatter(
            modeling.oof_predictions[observed_col],
            modeling.oof_predictions[predicted_col],
            s=14,
            alpha=0.45,
            color="#7B2CBF",
        )
        ax.set_xlabel(f"Observed {target}")
        ax.set_ylabel("Group-disjoint OOF prediction")
        ax.set_title(
            f"Best OOF model: {model}\nSpearman rho={float(best['spearman_rho']):.3f}"
        )
        fig.tight_layout()
        path = plot_dir / "06_best_oof_model.png"
        fig.savefig(path, dpi=180)
        plt.close(fig)
        written.append(path)
    return written


def file_sha256(path: Path) -> str:
    digest = sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def json_ready(value: object) -> object:
    if isinstance(value, Path):
        return str(value)
    if isinstance(value, (np.integer,)):
        return int(value)
    if isinstance(value, (np.floating,)):
        return None if not np.isfinite(value) else float(value)
    if isinstance(value, np.bool_):
        return bool(value)
    if isinstance(value, dict):
        return {str(key): json_ready(item) for key, item in value.items()}
    if isinstance(value, (list, tuple)):
        return [json_ready(item) for item in value]
    return value


def write_summary_markdown(
    path: Path,
    join_result: JoinResult,
    table: pd.DataFrame,
    eligible: pd.Series,
    category_column: str,
    category_summary: pd.DataFrame,
    correlations: pd.DataFrame,
    modeling: ModelingResult | None,
    warnings: list[str],
) -> None:
    lines = [
        "# FACS-guided 5'UTR retraining summary",
        "",
        f"- Script version: `{SCRIPT_VERSION}`",
        f"- Selected-library rows: {len(table) + len(join_result.unmatched_library):,}",
        f"- Sequence-matched rows: {len(table):,} ({join_result.match_rate:.2%})",
        f"- Analysis-eligible rows: {int(eligible.sum()):,}",
        f"- Exact `{category_column}` categories: {table[category_column].nunique(dropna=False)}",
        "",
        "## Category distribution",
        "",
        "```text",
        category_summary[["category", "selected_n", "eligible_n", "high15_top15_n", "high15_top15_rate"]]
        .to_string(index=False),
        "```",
        "",
        "## Original scores versus FACS endpoints",
        "",
    ]
    if len(correlations):
        top = correlations.sort_values("spearman_rho", key=lambda x: x.abs(), ascending=False).head(10)
        lines.extend(["```text", top[["old_score", "facs_target", "n", "spearman_rho", "mean_absolute_percentile_gap"]].to_string(index=False), "```"])
    else:
        lines.append("No usable original score columns were found.")
    lines.extend(["", "## Leakage-aware retraining", ""])
    if modeling is not None and len(modeling.metrics):
        lines.extend(
            [
                "Reported performance below is from gene/sequence-cluster-disjoint out-of-fold predictions.",
                "",
                "```text",
                modeling.metrics[["target", "model", "n", "n_features", "spearman_rho", "mae", "r2", "top15_roc_auc"]].to_string(index=False),
                "```",
            ]
        )
    else:
        lines.append("Modeling was skipped.")
    if warnings:
        lines.extend(["", "## Warnings", "", *[f"- {warning}" for warning in warnings]])
    lines.extend(
        [
            "",
            "## Interpretation boundary",
            "",
            "The experiment is one pooled Sort-seq run. Correlations and OOF performance quantify "
            "internal predictive consistency within this selected 2,000-UTR design space; they are "
            "not biological-replicate confidence intervals and do not establish translation-rate causality.",
            "Final joblib models are fit on all eligible rows only after OOF evaluation.",
        ]
    )
    path.write_text("\n".join(lines) + "\n", encoding="utf-8")


def run_analysis(args: argparse.Namespace) -> dict[str, object]:
    args.outdir.mkdir(parents=True, exist_ok=True)
    table_dir = args.outdir / "tables"
    table_dir.mkdir(parents=True, exist_ok=True)
    targets = [item.strip() for item in args.targets.split(",") if item.strip()]
    extra_scores = [item.strip() for item in args.old_score_columns.split(",") if item.strip()]
    fractions = [float(item) for item in args.top_fractions.split(",") if item.strip()]
    if not targets:
        raise ValueError("At least one FACS target is required")
    if any(not 0 < value < 1 for value in fractions):
        raise ValueError("All --top-fractions values must be between 0 and 1")

    library = read_table(args.selected_library)
    sortseq = read_table(args.sortseq_results)
    join_result = join_library_to_sortseq(
        library,
        sortseq,
        library_sequence_column=args.library_sequence_column,
        sortseq_sequence_column=args.sortseq_sequence_column,
        min_match_rate=0.0,
    )
    join_result.unmatched_library.to_csv(table_dir / "unmatched_selected_library.csv", index=False)
    join_result.unmatched_sortseq.to_csv(table_dir / "unmatched_sortseq_variants.csv", index=False)
    if join_result.match_rate < args.min_match_rate:
        raise ValueError(
            f"Sequence match rate {join_result.match_rate:.3%} is below --min-match-rate "
            f"{args.min_match_rate:.3%}. Inspect {table_dir / 'unmatched_selected_library.csv'} "
            "and unmatched_sortseq_variants.csv."
        )
    table, eligible, resolved_targets = add_analysis_annotations(
        join_result.joined,
        targets=targets,
        eligibility_column=args.eligibility_column,
        min_total_bin_count=args.min_total_bin_count,
    )
    if PRIMARY_TARGET not in resolved_targets:
        raise ValueError(
            f"The primary endpoint {PRIMARY_TARGET!r} is required for category hit enrichment"
        )
    category_column = resolve_selection_column(table, args.category_column)
    old_scores = choose_old_score_columns(table, extra_scores)
    correlations, category_correlations = score_correlations(
        table,
        eligible,
        old_scores,
        resolved_targets,
        category_column,
    )
    overlaps = topn_overlap(table, eligible, old_scores, resolved_targets)
    category_summary, category_tiers = category_tables(
        table, eligible, category_column, resolved_targets
    )
    hit_enrichment = category_hit_enrichment(
        table,
        eligible,
        category_column,
        resolved_targets[PRIMARY_TARGET],
        fractions,
    )

    modeling: ModelingResult | None = None
    if not args.skip_modeling:
        modeling = run_oof_models(
            table,
            eligible,
            resolved_targets,
            old_scores,
            sequence_column="_sequence_key",
            category_column=category_column,
            outdir=args.outdir,
            kmer_max=args.kmer_max,
            cv_folds=args.cv_folds,
            n_estimators=args.n_estimators,
            min_samples_leaf=args.min_samples_leaf,
            n_jobs=args.n_jobs,
            random_seed=args.random_seed,
            save_final_models=not args.no_final_models,
        )

    table.to_csv(table_dir / "facs_joined_selected_library.csv", index=False)
    category_summary.to_csv(table_dir / "category_distribution_and_scores.csv", index=False)
    category_tiers.to_csv(table_dir / "category_by_facs_high15_tier.csv", index=False)
    hit_enrichment.to_csv(table_dir / "category_top_hit_enrichment.csv", index=False)
    correlations.to_csv(table_dir / "original_score_vs_facs_correlations.csv", index=False)
    category_correlations.to_csv(
        table_dir / "original_score_vs_facs_correlations_by_category.csv", index=False
    )
    overlaps.to_csv(table_dir / "original_score_vs_facs_topn_overlap.csv", index=False)
    if modeling is not None:
        modeling.oof_predictions.to_csv(table_dir / "facs_model_oof_predictions.csv", index=False)
        modeling.metrics.to_csv(table_dir / "facs_model_oof_metrics.csv", index=False)
        modeling.category_metrics.to_csv(
            table_dir / "facs_model_oof_metrics_by_category.csv", index=False
        )

    warnings: list[str] = []
    category_count = int(table[category_column].nunique(dropna=False))
    if args.expected_category_count and category_count != args.expected_category_count:
        warnings.append(
            f"Observed {category_count} exact library groups, not the expected "
            f"{args.expected_category_count}. The output preserves the actual library_group labels."
        )
    if int(eligible.sum()) < 1000:
        warnings.append(
            f"Only {int(eligible.sum())} rows passed the experimental read-support filter; "
            "inspect category-specific eligibility before interpreting category hit rates."
        )
    if not old_scores:
        warnings.append("No usable original prediction-score columns were found in the selected library.")

    plot_paths: list[Path] = []
    if not args.no_plots:
        plot_paths = write_plots(
            args.outdir,
            table,
            eligible,
            category_column,
            resolved_targets,
            correlations,
            hit_enrichment,
            modeling,
        )

    write_summary_markdown(
        args.outdir / "FACS_RETRAINING_SUMMARY.md",
        join_result,
        table,
        eligible,
        category_column,
        category_summary,
        correlations,
        modeling,
        warnings,
    )

    manifest = {
        "script_version": SCRIPT_VERSION,
        "selected_library": str(args.selected_library.resolve()),
        "selected_library_sha256": file_sha256(args.selected_library),
        "sortseq_results": str(args.sortseq_results.resolve()),
        "sortseq_results_sha256": file_sha256(args.sortseq_results),
        "library_sequence_column": join_result.library_sequence_column,
        "sortseq_sequence_column": join_result.sortseq_sequence_column,
        "match_rate": join_result.match_rate,
        "matched_n": len(table),
        "eligible_n": int(eligible.sum()),
        "category_column": category_column,
        "category_count": category_count,
        "category_counts": table[category_column].fillna("NA").value_counts().to_dict(),
        "targets": resolved_targets,
        "old_score_columns": old_scores,
        "min_total_bin_count": args.min_total_bin_count,
        "eligibility_column": args.eligibility_column,
        "modeling": modeling.group_summary if modeling is not None else None,
        "final_models": [str(path) for path in modeling.final_model_paths] if modeling else [],
        "plots": [str(path) for path in plot_paths],
        "warnings": warnings,
    }
    (args.outdir / "run_manifest.json").write_text(
        json.dumps(json_ready(manifest), indent=2, ensure_ascii=False) + "\n",
        encoding="utf-8",
    )
    print(f"[SAVED] {args.outdir / 'FACS_RETRAINING_SUMMARY.md'}")
    print(f"[SAVED] {table_dir / 'original_score_vs_facs_correlations.csv'}")
    if modeling is not None:
        print(f"[SAVED] {table_dir / 'facs_model_oof_metrics.csv'}")
    return manifest


def main(argv: list[str] | None = None) -> int:
    args = parse_args(argv)
    run_analysis(args)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
