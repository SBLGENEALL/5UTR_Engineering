import argparse
import subprocess
import sys
from pathlib import Path


PIPELINE_STEPS = {
    0: "scripts/00_check_environment.py",
    1: "scripts/01_prepare_reference_genome.py",
    2: "scripts/02_annotate_cds.py",
    3: "scripts/03_annotate_atlas_tss.py",
    4: "scripts/04_extract_5utr_sequences.py",
    5: "scripts/05_map_rnaseq.py",
    6: "scripts/06_quantify_rnaseq.py",
    7: "scripts/07_map_riboseq.py",
    8: "scripts/08_quantify_riboseq.py",
    9: "scripts/09_map_proteomics.py",
    10: "scripts/10_calculate_te.py",
    11: "scripts/11_normalize_te_labels.py",
    12: "scripts/12_extract_5utr_features.py",
    13: "scripts/13_remove_similar_sequences.py",
    14: "scripts/14_train_rank_model.py",
    15: "scripts/15_score_candidate_library.py",
    16: "scripts/16_select_top_candidates.py",
}


LEGACY_NO_CONFIG_STEPS = {13, 14, 16}


def run_step(step_number: int, config_path: str):
    script = PIPELINE_STEPS[step_number]
    script_path = Path(script)

    if not script_path.exists():
        raise FileNotFoundError(f"Missing script: {script}")

    print(f"\n[STEP {step_number:02d}] Running {script}")

    if step_number in LEGACY_NO_CONFIG_STEPS:
        cmd = [sys.executable, str(script_path)]
    else:
        cmd = [sys.executable, str(script_path), "--config", config_path]

    result = subprocess.run(cmd, check=False)

    if result.returncode != 0:
        raise RuntimeError(f"Step {step_number:02d} failed: {script}")


def main():
    parser = argparse.ArgumentParser(
        description="Run selected steps of the 5UTR engineering pipeline."
    )
    parser.add_argument("--start", type=int, required=True, help="Start step number")
    parser.add_argument("--end", type=int, required=True, help="End step number")
    parser.add_argument(
        "--config",
        default="config/config.yaml",
        help="Path to config YAML file",
    )

    args = parser.parse_args()

    if args.start < 0 or args.end > 16 or args.start > args.end:
        raise ValueError("Invalid step range. Use --start 0 --end 16.")

    for step in range(args.start, args.end + 1):
        run_step(step, args.config)

    print("\nPipeline finished successfully.")


if __name__ == "__main__":
    main()