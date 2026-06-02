# Archive policy — final numbered release

The final team-facing scripts are strictly numbered `00`–`10` in `01_pipeline/scripts/`.

Anything not used by the final numbered workflow should be stored under `99_archive/` rather than deleted.

## Main scripts to keep

```text
00_check_inputs.py
01_build_utr_database.py
02_tss_correction.py
03_map_rna_ribo_public_te.py
04_preprocess_heffner_proteomics.py
05_integrate_proteomics_multiomics.py
06_plot_multiomics_distributions.py
07_heavy_rnafold_kmer6_automl.py
08_jaccard_sequence_cluster_qc.py
09_cluster_aware_classification_benchmark.py
10_select_2000_cluster_diverse_library.py
common.py
run_00_full_final_pipeline.py
run_01_annotation_tss_publicTE.py
run_02_proteomics_multiomics.py
run_03_model_jaccard_select2000.py
```

## Archived examples

Legacy names such as `04b_*`, `05b_*`, `06b_*`, `09c_*`, `18_*`, `21_*`, `22_*`, `23_*`, and `24_*` were kept only under:

```text
99_archive/legacy_scripts_pre_numbering/
```

These are not needed for the final team workflow.
