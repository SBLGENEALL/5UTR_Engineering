from pathlib import Path
import argparse
import shutil

ROOT = Path.cwd()
ARCHIVE = ROOT / "99_archive" / "auto_archived_nonfinal"
KEEP_SCRIPTS = {
    "00_check_inputs.py",
    "01_build_utr_database.py",
    "02_tss_correction.py",
    "03_map_rna_ribo_public_te.py",
    "04_preprocess_heffner_proteomics.py",
    "05_integrate_proteomics_multiomics.py",
    "06_plot_multiomics_distributions.py",
    "07_heavy_rnafold_kmer6_automl.py",
    "08_jaccard_sequence_cluster_qc.py",
    "09_cluster_aware_classification_benchmark.py",
    "10_select_2000_cluster_diverse_library.py",
    "common.py",
    "run_00_full_final_pipeline.py",
    "run_01_annotation_tss_publicTE.py",
    "run_02_proteomics_multiomics.py",
    "run_03_model_jaccard_select2000.py",
}

def archive_file(path, dry_run):
    rel = path.relative_to(ROOT)
    dest = ARCHIVE / rel
    print(f"ARCHIVE {rel} -> {dest.relative_to(ROOT)}")
    if not dry_run:
        dest.parent.mkdir(parents=True, exist_ok=True)
        shutil.move(str(path), str(dest))

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--apply', action='store_true')
    ap.add_argument('--dry-run', action='store_true')
    args = ap.parse_args()
    dry = not args.apply

    script_dir = ROOT / '01_pipeline' / 'scripts'
    if script_dir.exists():
        for p in script_dir.glob('*.py'):
            if p.name not in KEEP_SCRIPTS:
                archive_file(p, dry)

    # Archive root-level old run bats except final main.
    for p in ROOT.glob('RUN_*.bat'):
        if p.name != 'RUN_FINAL_MAIN.bat':
            archive_file(p, dry)

    print('Dry run complete.' if dry else 'Archive complete.')

if __name__ == '__main__':
    main()
