from pathlib import Path
import subprocess, sys
steps=["00_check_inputs.py","01_build_utr_database.py","02_tss_correction.py","03_map_rna_ribo_public_te.py"]
script_dir=Path("01_pipeline")/"scripts"
for step in steps:
    print("\n"+"="*80); print("RUN", step); print("="*80)
    r=subprocess.run([sys.executable, str(script_dir/step)])
    if r.returncode: raise SystemExit(f"FAILED: {step}")
print("DONE: annotation + TSS correction + public TE labeling")
