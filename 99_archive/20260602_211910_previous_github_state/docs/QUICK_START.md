# QUICK START

## 1. Prepare raw data

Populate the 00_raw_data directory using the reference datasets described in the script catalog workbook.

## 2. Run the pipeline

Recommended execution order:

RUN_01_base_publicTE.bat
RUN_02_proteomics_multiomics.bat
RUN_03_automl_quick.bat
RUN_04_plot_proteomics.bat
RUN_05_final_libraries.bat

Or execute:

RUN_ALL.bat

## 3. Key output

Final library:
selected_2000_50_100bp_cluster_diverse_evidence_balanced_library.csv

## 4. Workstation deployment

For offline workstation environments, manually copy all required datasets into 00_raw_data before execution.