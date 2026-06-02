# Required raw inputs

This final numbered release assumes all raw data are supplied. TSS atlas is required.

```text
00_raw_data/01_ncbi_genome_annotation/
  GCF_003668045.1_CriGri-PICR_genomic.fna.gz
  GCF_003668045.1_CriGri-PICR_genomic.gff.gz

00_raw_data/02_tss_atlas_picr3/
  GSE159044_eTSS_NCBI_picr.bed.gz
  GSE159044_eTSS_NCBI_picr.meta.tsv.gz

00_raw_data/03_gse79512_rna_ribo/
  GSE79512_RNASeq_rawCount.txt.gz
  GSE79512_RiboSeq_rawCount.txt.gz

00_raw_data/05_cho_proteomics/
  Heffner_minimal.tsv
  ncbi_gene_mapping/gene2accession.gz
  ncbi_gene_mapping/gene_info.gz
```

`data_assets/Heffner_minimal.tsv.gz` is included. `RUN_FINAL_MAIN.sh` automatically unzips it to `00_raw_data/05_cho_proteomics/Heffner_minimal.tsv` if missing.
