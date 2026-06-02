#!/usr/bin/env bash
set -euo pipefail

# CHO5UTR final numbered pipeline.
# Run from repository root.
# Required raw inputs are listed in docs/INPUTS_REQUIRED_TSS_RAWDATA.md.

if [ -f data_assets/Heffner_minimal.tsv.gz ] && [ ! -f 00_raw_data/05_cho_proteomics/Heffner_minimal.tsv ]; then
  mkdir -p 00_raw_data/05_cho_proteomics
  gunzip -c data_assets/Heffner_minimal.tsv.gz > 00_raw_data/05_cho_proteomics/Heffner_minimal.tsv
fi

python 01_pipeline/scripts/run_00_full_final_pipeline.py
