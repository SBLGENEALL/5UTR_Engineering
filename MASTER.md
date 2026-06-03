# 5UTR Engineering MASTER

Last Updated: 2026-06-04
Source of Truth: GitHub `main` branch + this `MASTER.md`
Status: ACTIVE

---

## 1. Project Purpose

This repository is the official, team-shareable CHO 5' UTR engineering pipeline for discovering candidate 5' UTR sequences using CHO genome/GFF annotation, PICR3 TSS atlas correction, public RNA/Ribo TE labels, proteomics integration, RNAfold/k-mer modeling, sequence-cluster QC, and final library selection.

The current official release is designed to start from raw data only and produce a final 2,000-member cluster-diverse 50-100 bp 5' UTR library.

---

## 2. Official Baseline

Repository: `SBLGENEALL/5UTR_Engineering`
Branch: `main`
Release concept: Final numbered release
Baseline package/reference: `CHO5UTR_FINAL_RELEASE_GitHub_ready.zip`

Official run command:

```bash
conda activate /home/MCET03/conda_envs/utr_env
bash RUN_FINAL_MAIN.sh
```

Alternative direct command:

```bash
python 01_pipeline/scripts/run_00_full_final_pipeline.py
```

---

## 3. Required Raw Inputs

The pipeline starts from raw/reference inputs only:

1. NCBI CHO genome FASTA + GFF
2. PICR3 TSS atlas BED + metadata
3. GSE79512 RNA-seq raw counts
4. GSE79512 Ribo-seq raw counts
5. Heffner CHO proteomics minimal TSV
6. NCBI `gene2accession.gz` and `gene_info.gz`

---

## 4. Current Pipeline Flow

```text
CHO genome/GFF annotation
→ TSS-atlas correction
→ RNA/Ribo public TE labeling
→ Heffner proteomics mapping
→ multi-omics label generation
→ RNAfold/k-mer6/tree2000 heavy modeling
→ Jaccard sequence cluster QC
→ cluster-aware classification benchmark
→ final 2,000 cluster-diverse 50-100 bp 5' UTR library
```

---

## 5. Numbered Main Scripts

| Step | Script | Purpose |
|---:|---|---|
| 00 | `00_check_inputs.py` | Validate required raw inputs. TSS atlas is required in this final release. |
| 01 | `01_build_utr_database.py` | Build annotation-derived CHO 5' UTR database from NCBI FASTA/GFF. |
| 02 | `02_tss_correction.py` | Correct/support UTRs using PICR3 TSS atlas. |
| 03 | `03_map_rna_ribo_public_te.py` | Map RNA/Ribo counts and create robust public TE labels. |
| 04 | `04_preprocess_heffner_proteomics.py` | Process Heffner minimal TSV and map proteins to genes. |
| 05 | `05_integrate_proteomics_multiomics.py` | Add protein abundance/residual labels to UTR rows. |
| 06 | `06_plot_multiomics_distributions.py` | Generate TE/proteomics/multiomics distribution QC plots. |
| 07 | `07_heavy_rnafold_kmer6_automl.py` | Run heavy RNAfold/k-mer6/tree model benchmark. |
| 08 | `08_jaccard_sequence_cluster_qc.py` | Remove/cluster exact and near-duplicate sequences. |
| 09 | `09_cluster_aware_classification_benchmark.py` | Evaluate classification with gene/sequence-cluster-aware splits. |
| 10 | `10_select_2000_cluster_diverse_library.py` | Select final 2,000 evidence-balanced, cluster-diverse candidates. |

---

## 6. Final Outputs

```text
07_library_design/tables/selected_2000_50_100bp_cluster_diverse_evidence_balanced_library.csv
07_library_design/fasta/selected_2000_50_100bp_cluster_diverse_evidence_balanced_library.fasta
```

---

## 7. Active Development Notes

Current active script/fix from recent workspace discussions:

`03_map_rna_ribo_robust_public_te.py` / step 03 public TE mapping logic.

Auto sample assignment rule when sample groups are set to `auto`:

```text
rna_day3  -> s01-s03
ribo_day3 -> s04-s06
rna_day6  -> s07-s09
ribo_day6 -> s10-s12
```

This rule should be preserved or reconciled with the current step-03 script naming in the repository.

---

## 8. Project Management Rules

- GitHub `main` is the official project state.
- `MASTER.md` is the lightweight project consensus document.
- Workspace chats are disposable working areas.
- When a workspace becomes slow, create the next workspace: `5UTR_WORKSPACE_02`, `5UTR_WORKSPACE_03`, etc.
- At the end of a workspace, update only validated decisions into `MASTER.md`, `CHANGELOG.md`, and `NEXT_ACTIONS.md`.
- Failed attempts and long discussions should not be copied into `MASTER.md` unless they affect current project decisions.

---

## 9. New Chat Start Template

Use this from PC or mobile:

```text
5UTR 프로젝트.
GitHub main의 MASTER.md 기준으로 진행해줘.
오늘 작업:
- ...
```

---

## 10. Current Next Focus

- Confirm step-03 auto sample assignment against the current repository script name and implementation.
- Keep the final numbered pipeline stable and team-shareable.
- Use `CHANGELOG.md` for dated changes and `NEXT_ACTIONS.md` for the next practical tasks.
