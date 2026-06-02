# Workstation Deployment

## Environment assumptions

- Offline workstation
- Local Python environment
- ViennaRNA RNAfold installed
- Large-memory execution available

## Data placement

Copy all datasets into:

00_raw_data/

Expected categories include:

01_ncbi_genome_annotation
02_tss_atlas_picr3
03_gse79512_rna_ribo
04_reference_controls
05_cho_proteomics
06_naturecomm_2021

## Recommended execution

1. Build UTR database
2. TSS correction
3. RNA/Ribo mapping
4. Proteomics integration
5. Feature extraction
6. Jaccard clustering
7. Cluster-aware benchmark
8. Final library selection

## Final output

selected_2000_50_100bp_cluster_diverse_evidence_balanced_library.csv