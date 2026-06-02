# 5UTR Engineering Latest Pipeline Specification

This document fixes the latest intended workflow for the team-distribution repository.

Important rule:

- Do not overwrite the final repository with old scripts from the home-PC Archive folder.
- The home-PC Archive folder is reference-only.
- The final repository should keep one script per experimental/computational step.
- Each script must be runnable independently and also through `run_all.py` / `run_from_step.py`.

---

## Final workflow

```text
CHO reference genome/GFF
→ CDS annotation
→ Atlas TSS annotation
→ 5UTR extraction
→ RNA-seq count mapping
→ Ribo-seq count mapping
→ public TE construction
→ proteomics mapping
→ protein/RNA residual and percentile rank labels
→ 5UTR feature extraction
→ Jaccard/k-mer similarity filtering
→ rank-model training/benchmark
→ candidate library scoring
→ top 2000 cluster-diverse candidate selection
```

---

## Canonical scripts

| Step | Script | Status | Main role |
|---:|---|---|---|
| 00 | `00_check_environment.py` | placeholder | Check Python packages and external tools |
| 01 | `01_prepare_reference_genome.py` | to implement | Validate/register CriGri-PICR genome FASTA and GFF/GTF |
| 02 | `02_annotate_cds.py` | to implement | Parse CDS/start-codon/gene annotation from GFF/GTF |
| 03 | `03_annotate_atlas_tss.py` | to implement | Map GSE159044 atlas TSS BED to NCBI/PICR gene annotation |
| 04 | `04_extract_5utr_sequences.py` | to implement | Extract TSS-corrected 5UTR sequences from genome |
| 05 | `05_map_rnaseq.py` | to implement | Load/map public RNA-seq count table |
| 06 | `06_quantify_rnaseq.py` | to implement | Standardize RNA abundance per gene |
| 07 | `07_map_riboseq.py` | to implement | Load/map public Ribo-seq count table |
| 08 | `08_quantify_riboseq.py` | to implement | Standardize ribosome footprint abundance per gene |
| 09 | `09_map_proteomics.py` | to implement | Map Heffner 2020 CHO proteomics data to genes |
| 10 | `10_calculate_te.py` | implemented baseline | Calculate RNA/Ribo/protein-derived TE and residual metrics |
| 11 | `11_normalize_te_labels.py` | implemented baseline | Convert TE/residual metric into percentile/rank labels |
| 12 | `12_extract_5utr_features.py` | implemented baseline | Generate sequence-level 5UTR features |
| 13 | `13_remove_similar_sequences.py` | legacy implemented | Jaccard/k-mer similarity clustering/filtering |
| 14 | `14_train_rank_model.py` | legacy implemented | Cluster-aware rank-model benchmark/training |
| 15 | `15_score_candidate_library.py` | implemented baseline | Score candidate 5UTR library |
| 16 | `16_select_top_candidates.py` | legacy implemented | Select final top 2000 diverse candidates |

---

## Final data sources

Expected raw/reference data sources:

```text
data/reference/GCF_003668045.1_CriGri-PICR_genomic.fna.gz
data/reference/GCF_003668045.1_CriGri-PICR_genomic.gff.gz

data/raw/tss_atlas/GSE159044_eTSS_NCBI_picr.bed.gz
data/raw/tss_atlas/GSE159044_eTSS_NCBI_picr.meta.tsv.gz

data/raw/rna_ribo/GSE79512_RiboSeq_rawCount.txt.gz

data/raw/proteomics/Heffner_2020_CHO_hamster_proteomics_supp_table1.xlsx
```

Processed canonical files should be:

```text
data/processed/cds_annotation.csv
data/processed/tss_annotation.csv
data/processed/utr_sequences.csv
data/processed/rna_abundance.csv
data/processed/ribo_abundance.csv
data/processed/protein_abundance.csv
data/processed/te_metrics.csv
data/processed/te_rank_labels.csv
data/processed/feature_matrix.csv
```

Final model/design outputs:

```text
results/14_train_rank_model/model.pkl
results/15_score_candidate_library/candidate_scores.csv
results/16_select_top_candidates/top_2000_candidates.csv
```

---

## Archive policy

Old scripts from previous experiments should not be deleted immediately. Store them under:

```text
archive_do_not_run/legacy_home_pc_archive/
archive_do_not_run/source_zip/
```

Only the numbered scripts in `scripts/` should be considered runnable team-facing scripts.

---

## Current implementation priority

1. Restore/implement Steps 01-04 from the latest conversation logic.
2. Restore/implement Steps 05-09 using the public RNA/Ribo and proteomics sources.
3. Replace baseline Steps 10-12/15 with final project-specific implementations if needed.
4. Keep Step 13-16 aligned with the final Jaccard/rank-model/top-2000 selection logic.
