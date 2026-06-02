# CHO 5′ UTR Engineering Pipeline — Final Numbered Release

This repository is the clean, team-shareable pipeline for CHO 5′ UTR candidate discovery.

The final release is designed to start from raw data only:

1. NCBI CHO genome FASTA + GFF
2. PICR3 TSS atlas BED + metadata
3. GSE79512 RNA-seq raw counts
4. GSE79512 Ribo-seq raw counts
5. Heffner CHO proteomics minimal TSV
6. NCBI `gene2accession.gz` and `gene_info.gz`

The pipeline performs:

```text
CHO genome/GFF annotation
→ TSS-atlas correction
→ RNA/Ribo public TE labeling
→ Heffner proteomics mapping
→ multi-omics label generation
→ RNAfold/k-mer6/tree2000 heavy modeling
→ Jaccard sequence cluster QC
→ cluster-aware classification benchmark
→ final 2,000 cluster-diverse 50–100 bp 5′ UTR library
```

## Run everything

```bash
conda activate /home/MCET03/conda_envs/utr_env
bash RUN_FINAL_MAIN.sh
```

or directly:

```bash
python 01_pipeline/scripts/run_00_full_final_pipeline.py
```

## Numbered main scripts

| Step | Script | Purpose |
|---:|---|---|
| 00 | `00_check_inputs.py` | Validate required raw inputs. TSS atlas is required in this final release. |
| 01 | `01_build_utr_database.py` | Build annotation-derived CHO 5′ UTR database from NCBI FASTA/GFF. |
| 02 | `02_tss_correction.py` | Correct/support UTRs using PICR3 TSS atlas. |
| 03 | `03_map_rna_ribo_public_te.py` | Map RNA/Ribo counts and create robust public TE labels. |
| 04 | `04_preprocess_heffner_proteomics.py` | Process Heffner minimal TSV and map proteins to genes. |
| 05 | `05_integrate_proteomics_multiomics.py` | Add protein abundance/residual labels to UTR rows. |
| 06 | `06_plot_multiomics_distributions.py` | Generate TE/proteomics/multiomics distribution QC plots. |
| 07 | `07_heavy_rnafold_kmer6_automl.py` | Run heavy RNAfold/k-mer6/tree model benchmark. |
| 08 | `08_jaccard_sequence_cluster_qc.py` | Remove/cluster exact and near-duplicate sequences. |
| 09 | `09_cluster_aware_classification_benchmark.py` | Evaluate classification with gene/sequence-cluster-aware splits. |
| 10 | `10_select_2000_cluster_diverse_library.py` | Select final 2,000 evidence-balanced, cluster-diverse candidates. |

## Final outputs

```text
07_library_design/tables/selected_2000_50_100bp_cluster_diverse_evidence_balanced_library.csv
07_library_design/fasta/selected_2000_50_100bp_cluster_diverse_evidence_balanced_library.fasta
```

See `docs/` for runbook, script catalog, input manifest, and interpretation guide.
