# CHO 5′UTR Engineering Platform

This repository is the **v1.0 final release package** for the CHO 5′UTR engineering project.

The repository now uses the uploaded final release package as the canonical source:

```text
CHO5UTR_FINAL_RELEASE_GitHub_ready.zip
CHO5UTR_FINAL_script_catalog_runbook.xlsx
```

Do not use the previous ChatGPT reconstruction scripts as the reference workflow. The ZIP package is the authoritative team-distribution version.

---

## What this pipeline does

This platform supports reproducible CHO 5′UTR candidate discovery by combining:

1. TSS-corrected CHO 5′UTR database construction
2. RNA-seq and Ribo-seq mapping
3. Public TE / Ribo-TE calculation
4. Proteomics and multi-omics integration
5. RNAfold + k-mer feature extraction
6. Cluster-aware model benchmarking
7. Evidence-balanced, cluster-diverse final library selection

The final goal is to generate an experimentally usable 5′UTR candidate library, especially:

```text
selected_2000_50_100bp_cluster_diverse_evidence_balanced_library.csv
```

---

## Final release package

Download and unzip:

```text
CHO5UTR_FINAL_RELEASE_GitHub_ready.zip
```

The unzipped release contains the actual executable pipeline and final folder structure.

Use the Excel runbook for script-by-script operation details:

```text
CHO5UTR_FINAL_script_catalog_runbook.xlsx
```

---

## Canonical release structure

After unzipping the release package, the expected project structure is:

```text
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
```

### Key folders

```text
00_raw_data/          Raw input data location
01_pipeline/         Main executable scripts and run scripts
02_utr_database/     UTR database construction outputs
03_tss_correction/   TSS-corrected 5′UTR outputs
04_te_labeling/      TE / protein residual / evidence labeling
05_feature_extraction/ RNAfold, k-mer, and sequence feature outputs
06_modeling/         Model training and cluster-aware benchmark outputs
07_library_design/   Final selected 5′UTR library outputs
08_reports/          Reports, plots, summaries, and runbooks
99_archive/          Deprecated or historical scripts; do not run by default
```

---

## Quick start

Use the batch files included in the final ZIP package.

Recommended order:

```text
RUN_00_download.bat
RUN_01_base_publicTE.bat
RUN_02_proteomics_multiomics.bat
RUN_03_automl_quick.bat
RUN_04_plot_proteomics.bat
RUN_05_final_libraries.bat
```

Or run the full workflow:

```text
RUN_ALL.bat
```

For workstation/server environments without internet access, prepare `00_raw_data/` manually before running the pipeline.

---

## Key scripts

The most important final-release scripts are inside:

```text
01_pipeline/scripts/
```

Core scripts:

```text
01_build_utr_database.py
02_tss_correction.py
03_map_rna_ribo_robust_public_te.py
09_integrate_proteomics_multiomics.py
18_heavy_rnafold_kmer6_automl.py
22_jaccard_sequence_cluster_qc.py
23_cluster_aware_classification_benchmark.py
24_select_2000_cluster_diverse_library.py
```

The final library selection script is:

```text
01_pipeline/scripts/24_select_2000_cluster_diverse_library.py
```

---

## Final output

Primary final candidate library:

```text
07_library_design/selected_2000_50_100bp_cluster_diverse_evidence_balanced_library.csv
```

This file is designed to balance multiple evidence groups, including:

```text
publicTE_high_confidence
TE_model_classifier_supported
Protein_abundance_supported
Protein_residual_supported
Multi_omics_consensus_high
sequence_diverse_exploratory
Length_GC_uAUG_diversity
low_signal_negative_controls
fill_best_remaining_allow_cluster2
```

---

## Archive policy

Scripts in `99_archive/` are preserved for traceability but should not be used as the default workflow.

The current default workflow is the v1.0 final release ZIP package plus the script catalog Excel file.

---

## Version

Current release:

```text
v1.0_final_release
```
