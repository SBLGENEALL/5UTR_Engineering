import csv
import subprocess
import sys
import tempfile
import unittest
from pathlib import Path

import pandas as pd


REPO = Path(__file__).resolve().parents[1]
SCRIPT = REPO / "01_pipeline/scripts/10_select_2000_cluster_diverse_library.py"


def encoded_sequence(number):
    alphabet = "ACG"
    suffix = ""
    for _ in range(10):
        suffix = alphabet[number % 3] + suffix
        number //= 3
    return "AC" * 25 + suffix


class V14SelectionRefillTests(unittest.TestCase):
    def test_non_j_k_refill_reaches_2000_with_caps(self):
        with tempfile.TemporaryDirectory() as tmp:
            workdir = Path(tmp)
            input_csv = workdir / "candidates.csv"
            fields = [
                "utr_id",
                "gene_name",
                "utr5_sequence_tss_corrected",
                "seq_cluster_id",
                "robust_public_te_rank",
                "heavy_ensemble_score",
                "protein_abundance_rank",
                "protein_residual_rank",
                "multi_omics_utr_rank",
                "is_expressed_public",
            ]
            with input_csv.open("w", newline="", encoding="utf-8") as handle:
                writer = csv.DictWriter(handle, fieldnames=fields)
                writer.writeheader()
                for i in range(2200):
                    writer.writerow(
                        {
                            "utr_id": f"utr_{i}",
                            "gene_name": f"gene_{i // 3}",
                            "utr5_sequence_tss_corrected": encoded_sequence(i),
                            "seq_cluster_id": f"cluster_{i // 2}",
                            "robust_public_te_rank": 0.5 + i / 2000 if i < 800 else "",
                            "heavy_ensemble_score": 0.5 + i / 2400 if i < 1200 else "",
                            "protein_abundance_rank": 0.7 if 800 <= i < 1500 else "",
                            "protein_residual_rank": 0.7 if 900 <= i < 1550 else "",
                            "multi_omics_utr_rank": 0.7 if i < 650 else "",
                            "is_expressed_public": True,
                        }
                    )

            subprocess.run(
                [
                    sys.executable,
                    str(SCRIPT),
                    "--input",
                    str(input_csv),
                    "--n",
                    "2000",
                    "--max-per-cluster",
                    "1",
                    "--allow-cluster-fill",
                    "2",
                    "--max-per-gene",
                    "3",
                ],
                cwd=workdir,
                check=True,
                capture_output=True,
                text=True,
            )

            qc_path = workdir / "07_library_design/tables/v1.4_selection_policy_qc.csv"
            qc = pd.read_csv(qc_path).iloc[0]
            self.assertEqual(int(qc["selected_n"]), 2000)
            self.assertEqual(int(qc["shortage_n"]), 0)
            self.assertEqual(int(qc["J_fill_selected_n"]), 0)
            self.assertGreater(int(qc["K_count"]), 0)
            self.assertLessEqual(int(qc["max_per_gene"]), 3)
            self.assertLessEqual(int(qc["max_per_seq_cluster"]), 2)


if __name__ == "__main__":
    unittest.main()
