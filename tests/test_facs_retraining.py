import argparse
import importlib.util
from pathlib import Path
import sys
import tempfile
import unittest

import numpy as np
import pandas as pd


SCRIPT = (
    Path(__file__).resolve().parents[1]
    / "01_pipeline/scripts/12_retrain_with_facs_sortseq.py"
)
SPEC = importlib.util.spec_from_file_location("facs_retraining", SCRIPT)
MODULE = importlib.util.module_from_spec(SPEC)
sys.modules[SPEC.name] = MODULE
SPEC.loader.exec_module(MODULE)


def sequence_for(number: int, length: int = 60) -> str:
    state = number + 17
    letters = []
    for position in range(length):
        state = (1103515245 * state + 12345 + position) % (2**31)
        letters.append("ACGT"[(state >> 7) % 4])
    letters[number % (length - 3) : number % (length - 3) + 3] = list(
        "GCG" if number % 2 else "ATA"
    )
    return "".join(letters)


def fixture_tables(n: int = 64) -> tuple[pd.DataFrame, pd.DataFrame]:
    sequences = [sequence_for(index) for index in range(n)]
    gc = np.array(
        [(sequence.count("G") + sequence.count("C")) / len(sequence) for sequence in sequences]
    )
    motif = np.array([sequence.count("GCG") for sequence in sequences], dtype=float)
    old_score = pd.Series(gc + 0.03 * motif).rank(pct=True).to_numpy()
    high15 = np.clip(0.04 + 0.32 * old_score + 0.015 * (np.arange(n) % 5), 0, 0.95)
    expected = 1.5 + 4.0 * high15 + 0.05 * (np.arange(n) % 3)
    categories = [f"category_{index % 8 + 1}" for index in range(n)]
    library = pd.DataFrame(
        {
            "library_index": np.arange(1, n + 1),
            "library_group": categories,
            "utr_id": [f"utr_{index:03d}" for index in range(n)],
            "gene_id": [f"gene_{index // 2:03d}" for index in range(n)],
            "seq_cluster_id": [f"cluster_{index // 2:03d}" for index in range(n)],
            "utr5_sequence_tss_corrected": sequences,
            "heavy_ensemble_score": old_score,
            "model_support_score": old_score * 0.9 + 0.05,
            "robust_public_te_rank": old_score,
        }
    )
    sortseq = pd.DataFrame(
        {
            "variant_id": [*sequences, sequence_for(999)],
            "total_6bin_count": [*([500] * n), 600],
            "high_bin_raw_count": [*([80] * n), 100],
            "top15_read_support_pass": [*([True] * n), True],
            "pass_coverage": [*([True] * n), True],
            "high15_probability": [*high15, 0.15],
            "expected_bin_score": [*expected, 3.0],
        }
    )
    return library, sortseq


class FacsRetrainingTests(unittest.TestCase):
    def test_sequence_join_preserves_all_selected_variants(self):
        library, sortseq = fixture_tables()
        joined = MODULE.join_library_to_sortseq(library, sortseq, min_match_rate=0.99)
        self.assertEqual(len(joined.joined), len(library))
        self.assertEqual(len(joined.unmatched_library), 0)
        self.assertEqual(len(joined.unmatched_sortseq), 1)
        self.assertAlmostEqual(joined.match_rate, 1.0)
        self.assertEqual(joined.library_sequence_column, "utr5_sequence_tss_corrected")
        self.assertEqual(joined.sortseq_sequence_column, "variant_id")

    def test_category_distribution_and_hit_enrichment_use_actual_groups(self):
        library, sortseq = fixture_tables()
        joined = MODULE.join_library_to_sortseq(library, sortseq)
        table, eligible, targets = MODULE.add_analysis_annotations(
            joined.joined,
            targets=["high15_probability", "expected_bin_score"],
            eligibility_column="top15_read_support_pass",
            min_total_bin_count=200,
        )
        summary, tiers = MODULE.category_tables(
            table, eligible, "library_group", targets
        )
        enrichment = MODULE.category_hit_enrichment(
            table,
            eligible,
            "library_group",
            targets["high15_probability"],
            [0.15],
        )
        self.assertEqual(len(summary), 8)
        self.assertTrue((summary["selected_n"] == 8).all())
        self.assertEqual(tiers["category"].nunique(), 8)
        self.assertEqual(enrichment["category"].nunique(), 8)
        self.assertEqual(int(enrichment["observed_hits"].sum()), 10)

    def test_gene_and_sequence_cluster_are_disjoint_components(self):
        library, _ = fixture_tables()
        groups, summary = MODULE.build_leakage_groups(library)
        self.assertEqual(summary["n_components"], 32)
        self.assertEqual(summary["largest_component_n"], 2)
        for _, indices in library.groupby("gene_id").groups.items():
            self.assertEqual(groups.loc[list(indices)].nunique(), 1)
        for _, indices in library.groupby("seq_cluster_id").groups.items():
            self.assertEqual(groups.loc[list(indices)].nunique(), 1)

    def test_oof_models_emit_predictions_for_both_endpoints(self):
        library, sortseq = fixture_tables()
        joined = MODULE.join_library_to_sortseq(library, sortseq)
        table, eligible, targets = MODULE.add_analysis_annotations(
            joined.joined,
            targets=["high15_probability", "expected_bin_score"],
            eligibility_column="top15_read_support_pass",
            min_total_bin_count=200,
        )
        old_scores = MODULE.choose_old_score_columns(table)
        with tempfile.TemporaryDirectory() as tmpdir:
            result = MODULE.run_oof_models(
                table,
                eligible,
                targets,
                old_scores,
                sequence_column="_sequence_key",
                category_column="library_group",
                outdir=Path(tmpdir),
                kmer_max=2,
                cv_folds=4,
                n_estimators=20,
                min_samples_leaf=2,
                n_jobs=1,
                random_seed=7,
                save_final_models=False,
            )
        self.assertEqual(set(result.metrics["target"]), {"high15_probability", "expected_bin_score"})
        self.assertEqual(set(result.metrics["model"]), {"sequence_only", "prior_only", "hybrid"})
        self.assertEqual(len(result.oof_predictions), len(library))
        prediction_columns = [
            column for column in result.oof_predictions if column.startswith("oof__")
        ]
        self.assertEqual(len(prediction_columns), 6)
        self.assertFalse(result.oof_predictions[prediction_columns].isna().any().any())
        self.assertGreater(
            float(
                result.metrics.loc[
                    result.metrics["model"].eq("prior_only"), "spearman_rho"
                ].min()
            ),
            0.7,
        )

    def test_end_to_end_writes_required_audit_tables(self):
        library, sortseq = fixture_tables()
        with tempfile.TemporaryDirectory() as tmpdir:
            root = Path(tmpdir)
            library_path = root / "selected.csv"
            sortseq_path = root / "utr_results_full.tsv"
            outdir = root / "report"
            library.to_csv(library_path, index=False)
            sortseq.to_csv(sortseq_path, sep="\t", index=False)
            args = argparse.Namespace(
                selected_library=library_path,
                sortseq_results=sortseq_path,
                outdir=outdir,
                library_sequence_column=None,
                sortseq_sequence_column=None,
                category_column="library_group",
                eligibility_column="top15_read_support_pass",
                min_total_bin_count=200,
                min_match_rate=0.95,
                expected_category_count=8,
                targets="high15_probability,expected_bin_score",
                old_score_columns="",
                top_fractions="0.05,0.10,0.15,0.20",
                cv_folds=4,
                kmer_max=2,
                n_estimators=10,
                min_samples_leaf=2,
                n_jobs=1,
                random_seed=11,
                skip_modeling=False,
                no_final_models=True,
                no_plots=True,
            )
            manifest = MODULE.run_analysis(args)
            required = [
                outdir / "FACS_RETRAINING_SUMMARY.md",
                outdir / "run_manifest.json",
                outdir / "tables/category_distribution_and_scores.csv",
                outdir / "tables/category_top_hit_enrichment.csv",
                outdir / "tables/original_score_vs_facs_correlations.csv",
                outdir / "tables/facs_model_oof_metrics.csv",
                outdir / "tables/facs_model_oof_predictions.csv",
            ]
            self.assertTrue(all(path.exists() for path in required))
            self.assertEqual(manifest["category_count"], 8)
            self.assertEqual(manifest["matched_n"], 64)


if __name__ == "__main__":
    unittest.main()
