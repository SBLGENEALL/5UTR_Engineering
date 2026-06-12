from pathlib import Path
import subprocess
import sys

script_dir = Path("01_pipeline") / "scripts"
steps = [
    "00_check_inputs.py",
    "01_build_utr_database.py",
    "02_tss_correction.py",
    "03_map_rna_ribo_public_te.py",
    "04_preprocess_heffner_proteomics.py",
    "05_integrate_proteomics_multiomics.py",
    "06_plot_multiomics_distributions.py",
    "08_jaccard_sequence_cluster_qc.py",
    "07_heavy_rnafold_kmer6_automl.py",
    "09_cluster_aware_classification_benchmark.py",
    "10_select_2000_cluster_diverse_library.py",
]

for step in steps:
    print("\n" + "=" * 90)
    print("RUN", step)
    print("=" * 90)
    cmd = [sys.executable, str(script_dir / step)]
    # final library selection defaults
    if step == "08_jaccard_sequence_cluster_qc.py":
        cmd += ["--k", "6", "--jaccard-threshold", "0.85", "--containment-threshold", "0.90", "--cluster-scope", "all"]
    elif step == "07_heavy_rnafold_kmer6_automl.py":
        cmd += [
            "--split-modes", "random,gene_split,seq_cluster_split,gene_seq_cluster_split",
            "--train-cluster-representative-only",
        ]
    elif step == "09_cluster_aware_classification_benchmark.py":
        cmd += ["--length-min", "20", "--length-max", "500", "--kmax", "5", "--n-estimators", "1000"]
    elif step == "10_select_2000_cluster_diverse_library.py":
        cmd += [
            "--n", "2000",
            "--max-per-cluster", "1",
            "--allow-cluster-fill", "2",
            "--max-per-gene", "3",
        ]
    r = subprocess.run(cmd)
    if r.returncode != 0:
        raise SystemExit(f"FAILED: {step}")

print("\nDONE: CHO5UTR final numbered pipeline completed.")
print("Final CSV : 07_library_design/tables/selected_2000_50_100bp_cluster_diverse_evidence_balanced_library.csv")
print("Final FASTA: 07_library_design/fasta/selected_2000_50_100bp_cluster_diverse_evidence_balanced_library.fasta")
print("Run v1.4 QC audit separately after installing minimap2 or BLAST+:")
print("python 01_pipeline/scripts/run_04_v14_qc_audit.py --baseline-v1.3 <v1.3 selected CSV>")
