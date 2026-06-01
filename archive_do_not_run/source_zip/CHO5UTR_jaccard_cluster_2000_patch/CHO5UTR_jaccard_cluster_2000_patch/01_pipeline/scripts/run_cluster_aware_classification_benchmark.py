import subprocess, sys
cmd=[sys.executable,"01_pipeline/scripts/23_cluster_aware_classification_benchmark.py","--length-min","20","--length-max","500","--kmax","5","--n-estimators","1000"]
print("RUN"," ".join(cmd))
r=subprocess.run(cmd)
raise SystemExit(r.returncode)
