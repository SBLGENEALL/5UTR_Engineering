# Validation Plan

## PR1.1 Baseline Checks

Run after sequence clustering, heavy modeling, benchmark, and selection:

```bash
python 01_pipeline/scripts/09_cluster_aware_classification_benchmark.py --length-min 20 --length-max 500 --kmax 5 --n-estimators 1000
python 01_pipeline/scripts/10_select_2000_cluster_diverse_library.py --n 2000 --max-per-cluster 1 --allow-cluster-fill 2
```

Expected reports:

```text
06_modeling/tables/cluster_split_disjointness_check.csv
06_modeling/tables/final_library_gene_cluster_diversity_summary.txt
```

Required split checks:

```text
gene_split: gene overlap == 0
seq_cluster_split: seq_cluster_id overlap == 0
gene_seq_cluster_split: gene overlap == 0 and seq_cluster_id overlap == 0
pass_required_for_split: all True
```

## PR2 TSS And Expression Checks

Compile the changed scripts:

```bash
python -m py_compile 01_pipeline/scripts/02_tss_correction.py
python -m py_compile 01_pipeline/scripts/03_map_rna_ribo_public_te.py
python -m py_compile 01_pipeline/scripts/10_select_2000_cluster_diverse_library.py
```

Run the affected steps:

```bash
python 01_pipeline/scripts/02_tss_correction.py --max-extend 300
python 01_pipeline/scripts/03_map_rna_ribo_public_te.py
python 01_pipeline/scripts/08_jaccard_sequence_cluster_qc.py --k 6 --jaccard-threshold 0.85 --containment-threshold 0.90 --cluster-scope all
python 01_pipeline/scripts/07_heavy_rnafold_kmer6_automl.py --split-modes random,gene_split,seq_cluster_split,gene_seq_cluster_split --train-cluster-representative-only
python 01_pipeline/scripts/09_cluster_aware_classification_benchmark.py --length-min 20 --length-max 500 --kmax 5 --n-estimators 1000
python 01_pipeline/scripts/10_select_2000_cluster_diverse_library.py --n 2000 --max-per-cluster 1 --allow-cluster-fill 2
```

Expected PR2 reports:

```text
03_tss_correction/qc/tss_correction_summary.txt
04_te_labeling/qc/robust_public_te_mapping_summary.txt
07_library_design/qc/selected_2000_50_100bp_cluster_diverse_evidence_balanced_summary.txt
```

QC criteria:

```text
tss_correction_mode includes no_tss_match, trim_to_tss, extend_to_tss, or unchanged.
robust_public_te_rank is non-null only when is_expressed_public is True.
Non-expressed genes have expression_qc_reason explaining the missing or low-count group.
A/B/C/D/E library groups are selected from expressed evidence candidates.
F/G/H/J may use clean base candidates when evidence candidates are insufficient or not required by group design.
Selection summary reports candidate_pool_after_QC, evidence_candidate_pool_after_expression_TE_QC, and selection_source counts.
Final output filenames remain unchanged.
selected_n should remain 2000; if candidate_pool_after_QC drops below 2000, relax expression thresholds only after review.
```

## v1.4 PR1 QC Audit

Compile changed scripts:

```bash
python -m py_compile \
  01_pipeline/scripts/04_preprocess_heffner_proteomics.py \
  01_pipeline/scripts/05_integrate_proteomics_multiomics.py \
  01_pipeline/scripts/10_select_2000_cluster_diverse_library.py \
  01_pipeline/scripts/11_v14_qc_audit.py
```

Run affected stages:

```bash
python 01_pipeline/scripts/04_preprocess_heffner_proteomics.py
python 01_pipeline/scripts/05_integrate_proteomics_multiomics.py
python 01_pipeline/scripts/10_select_2000_cluster_diverse_library.py \
  --n 2000 --max-per-cluster 1 --allow-cluster-fill 2 --max-per-gene 3
python 01_pipeline/scripts/11_v14_qc_audit.py \
  --baseline-v1.3 07_library_design/tables/v1.3_selected_2000_library.csv
```

Required checks:

```text
accession-to-gene conversion success rate is reported.
gene-to-corrected-UTR mapping success rate is reported.
unmapped accession and unmapped gene tables are emitted.
selected_n == 2000.
shortage_n == 0.
J_fill_selected_n == 0.
selection refill uses ordered K1 A/B/E evidence, K2 C/D proteomics,
K3 classifier/model, and K4 F/G diversity pools.
F/G diversity and H negative-control groups remain available.
max_per_gene <= 3.
max_per_seq_cluster <= 2.
v1.4_selection_policy_qc.csv reports A/B/C/D/E/F/G/H/K counts and
protein-, classifier-, and multiomics-supported selected totals.
CHO mapping summary contains all selected library rows.
suspected non-CHO/hallucinated sequences are explicitly listed.
Excel, PPT, length/GC tables, charts, and v1.3 comparison are generated.
```
