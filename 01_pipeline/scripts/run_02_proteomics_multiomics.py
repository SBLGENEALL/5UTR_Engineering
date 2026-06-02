from pathlib import Path
import subprocess, sys
steps=["04_preprocess_heffner_proteomics.py","05_integrate_proteomics_multiomics.py","06_plot_multiomics_distributions.py"]
script_dir=Path("01_pipeline")/"scripts"
for step in steps:
    print("\n"+"="*80); print("RUN", step); print("="*80)
    r=subprocess.run([sys.executable, str(script_dir/step)])
    if r.returncode: raise SystemExit(f"FAILED: {step}")
print("DONE: proteomics + multiomics integration + distribution QC")
