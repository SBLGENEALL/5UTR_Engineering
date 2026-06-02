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


# Script catalog — final numbered release

| Step | Script | Input | Main output | Meaning |
|---:|---|---|---|---|
| 00 | `00_check_inputs.py` | raw data folders | terminal PASS/ERROR | Confirms all required raw inputs are present, including TSS atlas. |
| 01 | `01_build_utr_database.py` | NCBI FASTA/GFF | `02_utr_database/tables/` | Creates annotation-based CHO 5′ UTR candidates. |
| 02 | `02_tss_correction.py` | UTR candidates + PICR3 TSS atlas | `03_tss_correction/tables/tss_corrected_5utr_database.csv` | Supports/corrects 5′ UTRs using TSS signal. |
| 03 | `03_map_rna_ribo_public_te.py` | TSS-corrected UTRs + RNA/Ribo counts | `04_te_labeling/tables/tss_corrected_5utr_robust_public_te_labels.csv` | Creates robust public TE rank labels. |
| 04 | `04_preprocess_heffner_proteomics.py` | Heffner minimal TSV + NCBI gene mapping | `00_raw_data/05_cho_proteomics/*mapped_for_5utr.csv` | Produces protein abundance proxy mapped to gene IDs/symbols. |
| 05 | `05_integrate_proteomics_multiomics.py` | publicTE labels + mapped proteomics | `04_te_labeling/tables/tss_corrected_5utr_multiomics_labels.csv` | Adds protein abundance rank, protein residual rank, and multiomics rank. |
| 06 | `06_plot_multiomics_distributions.py` | multiomics labels | `06_modeling/plots/multiomics_distributions/` | QC for TE/protein/multiomics distributions and correlations. |
| 07 | `07_heavy_rnafold_kmer6_automl.py` | multiomics labels | heavy model benchmark tables/plots | Runs RNAfold/k-mer6/tree model search; classification is main interpretation. |
| 08 | `08_jaccard_sequence_cluster_qc.py` | multiomics labels | `04_te_labeling/tables/tss_corrected_5utr_with_seq_clusters.csv` | Clusters exact/near-duplicate UTRs by k-mer Jaccard/containment. |
| 09 | `09_cluster_aware_classification_benchmark.py` | sequence-clustered labels | `06_modeling/tables/cluster_aware_classification_benchmark.csv` | Tests high/low classification under gene/sequence-cluster splits. |
| 10 | `10_select_2000_cluster_diverse_library.py` | clustered labels + evidence scores | final CSV/FASTA in `07_library_design/` | Selects final 2,000 50–100 bp candidates with evidence + diversity balance. |

Legacy/basic 1000-candidate scripts and old patch runners are archived under `99_archive/legacy_scripts_pre_numbering/`.
