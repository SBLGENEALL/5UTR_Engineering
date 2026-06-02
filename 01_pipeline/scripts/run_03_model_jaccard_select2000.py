from pathlib import Path
import subprocess, sys
script_dir=Path("01_pipeline")/"scripts"
commands=[
    ["07_heavy_rnafold_kmer6_automl.py"],
    ["08_jaccard_sequence_cluster_qc.py","--k","6","--jaccard-threshold","0.85","--containment-threshold","0.90","--cluster-scope","all"],
    ["09_cluster_aware_classification_benchmark.py","--length-min","20","--length-max","500","--kmax","5","--n-estimators","1000"],
    ["10_select_2000_cluster_diverse_library.py","--n","2000","--max-per-cluster","1","--allow-cluster-fill","2"],
]
for cmd in commands:
    step=cmd[0]
    print("\n"+"="*80); print("RUN", " ".join(cmd)); print("="*80)
    r=subprocess.run([sys.executable, str(script_dir/step), *cmd[1:]])
    if r.returncode: raise SystemExit(f"FAILED: {step}")
print("DONE: heavy model + Jaccard + final 2000 selection")
