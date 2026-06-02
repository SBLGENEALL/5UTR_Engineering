import subprocess, sys
cmd=[sys.executable,"01_pipeline/scripts/22_jaccard_sequence_cluster_qc.py","--k","6","--jaccard-threshold","0.85","--containment-threshold","0.90","--cluster-scope","all"]
print("RUN"," ".join(cmd))
r=subprocess.run(cmd)
raise SystemExit(r.returncode)
