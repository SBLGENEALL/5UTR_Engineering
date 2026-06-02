# CHO5UTR Jaccard / cluster-aware validation / 2000 candidate library patch

This patch adds the next improvement stage:

1. Exact duplicate + 6-mer Jaccard/containment sequence clustering
2. Cluster-aware high/low classification benchmark
3. Evidence-balanced, cluster-diverse 2000-candidate 50-100 bp library selection

## Apply

From project root:

```bash
unzip -o CHO5UTR_jaccard_cluster_2000_patch.zip
```

## Recommended run order

```bash
python 01_pipeline/scripts/run_jaccard_cluster_qc.py
python 01_pipeline/scripts/run_cluster_aware_classification_benchmark.py
python 01_pipeline/scripts/run_select_2000_cluster_diverse_library.py
```

or all at once:

```bash
python 01_pipeline/scripts/run_jaccard_2000_full_pipeline.py
```

## Inputs

Preferred:

```text
04_te_labeling/tables/tss_corrected_5utr_multiomics_labels.csv
```

If unavailable, scripts fall back to:

```text
04_te_labeling/tables/tss_corrected_5utr_robust_public_te_labels.csv
```

## Outputs

### Jaccard clustering

```text
04_te_labeling/tables/tss_corrected_5utr_with_seq_clusters.csv
04_te_labeling/qc/jaccard_sequence_cluster_summary.csv
04_te_labeling/qc/jaccard_sequence_cluster_report.txt
```

### Cluster-aware benchmark

```text
06_modeling/tables/cluster_aware_classification_benchmark.csv
06_modeling/tables/cluster_aware_classification_summary.txt
06_modeling/plots/cluster_aware_classification_benchmark.png
```

### Final 2000 library

```text
07_library_design/tables/selected_2000_50_100bp_cluster_diverse_evidence_balanced_library.csv
07_library_design/fasta/selected_2000_50_100bp_cluster_diverse_evidence_balanced_library.fasta
07_library_design/qc/selected_2000_50_100bp_cluster_diverse_evidence_balanced_summary.txt
```

## Interpretation

Use cluster-aware classification as the more stringent benchmark:

- random split: optimistic reference only
- gene_split: gene-level generalization
- seq_cluster_split: sequence-similarity generalization
- gene_seq_cluster_split: strictest; no same gene or near-duplicate sequence leakage

Final 2000 library is not model top-2000. It is evidence-balanced:

- public TE high confidence
- TE/model-supported candidates
- protein abundance-supported candidates
- protein/RNA residual-supported candidates
- multiomics consensus candidates
- sequence-diverse exploratory candidates
- length/GC/uAUG diversity candidates
- negative controls
- optional reference controls if provided in 00_raw_data/04_reference_controls
