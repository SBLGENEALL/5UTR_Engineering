# 5UTR Engineering CHANGELOG

This file records only project-level decisions and validated changes. Long exploratory chat history should not be copied here.

---

## 2026-06-04

### Added

- Added `MASTER.md` as the lightweight project consensus document.
- Established GitHub `main` + `MASTER.md` as the source of truth for future ChatGPT workspaces.
- Added workspace rotation rule: use `5UTR_WORKSPACE_01`, `5UTR_WORKSPACE_02`, etc. and update `MASTER.md` before moving to the next workspace.

### Confirmed

- Repository is the final, team-shareable CHO 5' UTR candidate discovery pipeline.
- Final pipeline starts from raw/reference inputs only.
- Official run entry points:

```bash
bash RUN_FINAL_MAIN.sh
python 01_pipeline/scripts/run_00_full_final_pipeline.py
```

### Pipeline Consensus

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

### Active Development Note

- Step 03 public TE mapping / robust TE mapping should preserve or reconcile the recent auto sample assignment rule:

```text
rna_day3  -> s01-s03
ribo_day3 -> s04-s06
rna_day6  -> s07-s09
ribo_day6 -> s10-s12
```

### Final Outputs

```text
07_library_design/tables/selected_2000_50_100bp_cluster_diverse_evidence_balanced_library.csv
07_library_design/fasta/selected_2000_50_100bp_cluster_diverse_evidence_balanced_library.fasta
```
