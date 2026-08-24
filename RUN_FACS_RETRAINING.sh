#!/usr/bin/env bash
set -euo pipefail

ROOT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PYTHON_BIN="${PYTHON_BIN:-python}"

if [[ $# -lt 1 ]]; then
  cat <<'EOF'
Usage:
  bash RUN_FACS_RETRAINING.sh /absolute/path/to/results/sortseq/utr_results_full.tsv [extra options]

Example:
  FACS_N_JOBS=64 bash RUN_FACS_RETRAINING.sh \
    /data/user/MCET03/03_NGS/02_5UTR_sorting/results/sortseq/utr_results_full.tsv
EOF
  exit 2
fi

SORTSEQ_RESULTS="$1"
shift

export MPLCONFIGDIR="${MPLCONFIGDIR:-${TMPDIR:-/tmp}/cho5utr_matplotlib_cache}"
mkdir -p "$MPLCONFIGDIR"

"$PYTHON_BIN" "$ROOT_DIR/01_pipeline/scripts/12_retrain_with_facs_sortseq.py" \
  --selected-library "$ROOT_DIR/07_library_design/tables/selected_2000_50_100bp_cluster_diverse_evidence_balanced_library.csv" \
  --sortseq-results "$SORTSEQ_RESULTS" \
  --outdir "$ROOT_DIR/08_reports/facs_retraining" \
  --min-total-bin-count 200 \
  --eligibility-column top15_read_support_pass \
  --cv-folds 5 \
  --kmer-max 6 \
  --n-estimators 300 \
  --n-jobs "${FACS_N_JOBS:-32}" \
  "$@"
