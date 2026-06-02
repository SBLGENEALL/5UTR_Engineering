# Repository Structure

Canonical release layout:

00_raw_data/
01_pipeline/
02_utr_database/
03_tss_correction/
04_te_labeling/
05_feature_extraction/
06_modeling/
07_library_design/
08_reports/
99_archive/

## Important scripts

01_build_utr_database.py
02_tss_correction.py
03_map_rna_ribo_robust_public_te.py
09_integrate_proteomics_multiomics.py
18_heavy_rnafold_kmer6_automl.py
22_jaccard_sequence_cluster_qc.py
23_cluster_aware_classification_benchmark.py
24_select_2000_cluster_diverse_library.py

## Final deliverable

selected_2000_50_100bp_cluster_diverse_evidence_balanced_library.csv