import argparse
import subprocess
import sys
from pathlib import Path


DATA_BUILD_STEPS = [
    ["scripts/01_prepare_reference_genome.py"],
    ["scripts/02_annotate_cds.py"],
    ["scripts/03_annotate_atlas_tss.py"],
    ["scripts/04_extract_5utr_sequences.py"],
    ["scripts/05_map_rnaseq.py"],
    ["scripts/06_quantify_rnaseq.py"],
    ["scripts/07_map_riboseq.py"],
    ["scripts/08_quantify_riboseq.py"],
    ["scripts/09_map_proteomics.py"],
    ["scripts/10_calculate_te.py"],
    ["scripts/11_normalize_te_labels.py"],
]


def run_cmd(cmd):
    printable = " ".join(cmd)
    print(f"\n[RUN] {printable}", flush=True)
    result = subprocess.run(cmd, check=False)
    if result.returncode != 0:
        raise RuntimeError(f"Command failed: {printable}")


def ensure_required_files():
    required = [
        "data/reference/GCF_003668045.1_CriGri-PICR_genomic.fna.gz",
        "data/reference/GCF_003668045.1_CriGri-PICR_genomic.gff.gz",
        "data/raw/tss_atlas/GSE159044_eTSS_NCBI_picr.bed.gz",
        "data/raw/tss_atlas/GSE159044_eTSS_NCBI_picr.meta.tsv.gz",
        "data/raw/rna_ribo/GSE79512_RNASeq_rawCount.txt.gz",
        "data/raw/rna_ribo/GSE79512_RiboSeq_rawCount.txt.gz",
        "data/raw/proteomics/Heffner_2020_CHO_hamster_proteomics_ncbi_mapped_for_5utr.csv",
    ]
    missing = [p for p in required if not Path(p).exists()]
    if missing:
        raise FileNotFoundError("Missing required input files:\n" + "\n".join(missing))


def main():
    parser = argparse.ArgumentParser(description="Run final 5UTR workstation pipeline.")
    parser.add_argument("--mode", choices=["all", "data", "model"], default="all")
    parser.add_argument("--min-length", type=int, default=40)
    parser.add_argument("--max-length", type=int, default=200)
    parser.add_argument("--kmax", type=int, default=6)
    parser.add_argument("--use-rnafold", action="store_true", default=True)
    parser.add_argument("--skip-rnafold", action="store_true")
    parser.add_argument("--n-estimators", type=int, default=2000)
    parser.add_argument("--cluster-scope", default="all", choices=["all", "train_20_500", "selection_50_100"])
    args = parser.parse_args()

    ensure_required_files()

    if args.mode in ["all", "data"]:
        for step in DATA_BUILD_STEPS:
            run_cmd([sys.executable] + step)

    if args.mode in ["all", "model"]:
        feature_cmd = [
            sys.executable,
            "scripts/12_extract_5utr_features.py",
            "--min-length", str(args.min_length),
            "--max-length", str(args.max_length),
            "--kmax", str(args.kmax),
        ]
        if args.use_rnafold and not args.skip_rnafold:
            feature_cmd.append("--use-rnafold")
        run_cmd(feature_cmd)

        run_cmd([
            sys.executable,
            "scripts/13_remove_similar_sequences.py",
            "--input", "data/processed/feature_matrix.csv",
            "--cluster-scope", args.cluster_scope,
        ])

        run_cmd([
            sys.executable,
            "scripts/14_train_rank_model.py",
            "--input", "04_te_labeling/tables/tss_corrected_5utr_with_seq_clusters.csv",
            "--length-min", str(args.min_length),
            "--length-max", str(args.max_length),
            "--kmax", str(args.kmax),
            "--n-estimators", str(args.n_estimators),
        ])

    print("\n[DONE] Final workstation pipeline finished successfully.")


if __name__ == "__main__":
    main()
