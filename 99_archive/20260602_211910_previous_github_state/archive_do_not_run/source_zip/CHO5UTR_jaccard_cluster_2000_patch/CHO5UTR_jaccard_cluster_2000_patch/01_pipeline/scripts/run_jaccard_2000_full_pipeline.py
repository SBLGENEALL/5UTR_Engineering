import subprocess, sys
cmds=[
 [sys.executable,"01_pipeline/scripts/22_jaccard_sequence_cluster_qc.py","--k","6","--jaccard-threshold","0.85","--containment-threshold","0.90","--cluster-scope","all"],
 [sys.executable,"01_pipeline/scripts/23_cluster_aware_classification_benchmark.py","--length-min","20","--length-max","500","--kmax","5","--n-estimators","1000"],
 [sys.executable,"01_pipeline/scripts/24_select_2000_cluster_diverse_library.py","--n","2000","--max-per-cluster","1","--allow-cluster-fill","2"],
]
for cmd in cmds:
    print("\n"+"="*90)
    print("RUN"," ".join(cmd))
    print("="*90)
    r=subprocess.run(cmd)
    if r.returncode != 0:
        raise SystemExit(r.returncode)
