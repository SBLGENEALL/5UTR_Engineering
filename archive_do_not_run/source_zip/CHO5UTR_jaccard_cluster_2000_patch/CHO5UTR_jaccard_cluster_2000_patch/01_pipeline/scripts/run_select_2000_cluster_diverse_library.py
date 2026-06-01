import subprocess, sys
cmd=[sys.executable,"01_pipeline/scripts/24_select_2000_cluster_diverse_library.py","--n","2000","--max-per-cluster","1","--allow-cluster-fill","2"]
print("RUN"," ".join(cmd))
r=subprocess.run(cmd)
raise SystemExit(r.returncode)
