# Final pipeline order

## One-command execution

```bash
bash RUN_FINAL_MAIN.sh
```

## Step-by-step execution

```bash
python 01_pipeline/scripts/00_check_inputs.py
python 01_pipeline/scripts/01_build_utr_database.py
python 01_pipeline/scripts/02_tss_correction.py
python 01_pipeline/scripts/03_map_rna_ribo_public_te.py
python 01_pipeline/scripts/04_preprocess_heffner_proteomics.py
python 01_pipeline/scripts/05_integrate_proteomics_multiomics.py
python 01_pipeline/scripts/06_plot_multiomics_distributions.py
python 01_pipeline/scripts/07_heavy_rnafold_kmer6_automl.py
python 01_pipeline/scripts/08_jaccard_sequence_cluster_qc.py --k 6 --jaccard-threshold 0.85 --containment-threshold 0.90 --cluster-scope all
python 01_pipeline/scripts/09_cluster_aware_classification_benchmark.py --length-min 20 --length-max 500 --kmax 5 --n-estimators 1000
python 01_pipeline/scripts/10_select_2000_cluster_diverse_library.py --n 2000 --max-per-cluster 1 --allow-cluster-fill 2
```

## TSS atlas step

TSS correction is explicitly step 02:

```text
02_tss_correction.py
```

It uses:

```text
00_raw_data/02_tss_atlas_picr3/GSE159044_eTSS_NCBI_picr.bed.gz
00_raw_data/02_tss_atlas_picr3/GSE159044_eTSS_NCBI_picr.meta.tsv.gz
```

In this final release, the input checker treats TSS atlas as required because the team workflow will supply it.
