# CHANGELOG

## 2026-06-12 - v1.4 PR1 non-J refill correction

* Rebalanced A-H primary quotas and reduced H negative controls from 100 to 50.
* Added ordered K1-K4 non-J refill pools for A/B/E evidence, C/D proteomics,
  classifier/model support, and F/G diversity.
* Enforced final `max_per_gene <= 3` and `max_per_seq_cluster <= 2`.
* Added `v1.4_selection_policy_qc.csv` with shortage, group, and support totals.
* Kept uAUG0, 50-100 nt, forbidden-site, and unique-sequence filters unchanged.

## 2026-06-10 - v1.4 PR1 QC audit implementation

Branch:

* `v1.4-pr1-qc-audit`

Baseline:

* Created from `main`.
* Merged official validated `v1.3` into the PR branch.
* `main` itself was not modified or merged.

Changes:

* Added accession-to-gene and gene-to-corrected-UTR proteomics audit tables.
* Added unmapped accession and unmapped gene outputs.
* Removed generic `J_fill` production selection.
* Restricted final refill to A/B/C/E evidence, classifier, proteomics, and multiomics pools.
* Added CHO genome-origin validation using minimap2 or BLASTN.
* Added length/GC versus TE, heavy model, and protein QC.
* Added Excel, PPT, PNG chart, and v1.3-versus-v1.4 comparison outputs.

Policy:

* Evidence-only refill shortfalls fail by default.
* Raw sequence fallback is not silently relabeled as evidence.
* F/G diversity groups and H negative controls remain part of the library design.

## 2026-06-10 - v1.3 final library QC staged

Branch:

* `v1.3`

Status:

* PR3-2 workstation validation passed.
* v1.3 is the validated baseline to preserve before v1.4 development.
* No merge to `main` was performed.

Final QC:

```text
selected_n: 2000
requested_n: 2000
uaug_positive_n: 0
uaug0_policy_pass: True
n_unique_seq_clusters: 1937
max_per_seq_cluster: 2
cluster_cap: 2
cluster_cap_pass: True
gene_key: gene_name
n_unique_genes: 1896
max_per_gene: 3
gene_cap: 4
gene_cap_pass: True
mean_heavy_ensemble_score: 0.5811515
mean_robust_public_te_rank: 0.642380
```

Committed QC artifact:

```text
07_library_design/qc/v1.3_final_library_qc_summary.txt
```

Exclusions:

* Raw data.
* Intermediate tables and FASTA files.
* Caches and logs.
* Unvalidated generated outputs.

## 2026-06-05 - PR3-1 planned: uAUG source audit and uAUG=0 dry-run reporting

Branch:

* `improved-v1.2-pr3-selection-policy`

Baseline:

* `improved-v1.2`

Status:

* Reporting-only PR3-1 implementation.
* Production final library selection behavior is intentionally unchanged.
* Existing selected 2000 CSV/FASTA filenames remain unchanged.

Changed files:

```text
01_pipeline/scripts/10_select_2000_cluster_diverse_library.py
docs/CHANGELOG.md
```

Planned outputs:

```text
07_library_design/qc/uaug_source_by_group_summary.txt
07_library_design/tables/uaug_source_by_group_summary.csv
07_library_design/tables/uaug_positive_final_library_rows.csv
07_library_design/tables/uaug0_hard_filter_dry_run_summary.csv
07_library_design/tables/uaug0_hard_filter_quota_shortfall.csv
07_library_design/tables/uaug0_replacement_candidates.csv
```

Purpose:

* Identify which `library_group` and `selection_source` introduce uAUG-positive final library members.
* Quantify how many current final-library rows would be removed by a strict `uaug_count == 0` policy.
* Run a uAUG=0 hard-filter dry-run without overwriting the production selected CSV/FASTA.
* Report whether the 2,000-member library can be refilled from uAUG-free candidates while preserving quota, evidence, and sequence-cluster diversity constraints.

## 2026-06-04 - PR2 validated and promoted to improved-v1.2

Branch:

* `improved-v1.2`

Source branch:

* `improved-v1.1-pr2-tss-expression`

Status:

* Full raw-data validation passed on company workstation.
* `improved-v1.2` is now the validated working baseline for PR3 development.
* `main` remains the stable release branch and is not replaced yet.

Validation summary:

* Final library generated successfully.
* `selected_n = 2000`.
* CSV line count = 2001 including header.
* Length mean = 73.527.
* Length min/max = 50/100.
* `max_per_cluster_primary = 1`.
* `heavy_ensemble_score` and `heavy_ensemble_rank` are present in the final selected library.

TSS correction summary:

```text
trim_to_tss     22,722
extend_to_tss    9,799
unchanged         7,827
no_tss_match      2,112
total            42,460
```

TSS confidence summary:

```text
tss_supported_with_signal  40,348
no_tss_match                2,112
```

Model validation, strict `robust_public_te_rank` / `gene_seq_cluster_split` / `40_200` classification:

```text
ExtraTrees              ROC-AUC 0.7145
RandomForest            ROC-AUC 0.7092
HistGradientBoosting    ROC-AUC 0.6770
```

Reference PR1 baseline:

```text
RandomForest ROC-AUC 0.671
```

Leakage/disjointness validation:

* Strict `gene_seq_cluster_split` overlap/leakage was 0 for all four checked targets:
  * TE
  * multi-omics
  * protein abundance
  * protein residual

Final selection source counts:

```text
evidence_cand                 1238
fill_base_cand_clean_sequence  482
base_cand_diversity            100
base_cand_low_publicTE         100
fill_evidence_cand              68
base_cand_exploratory           12
```

uAUG audit observations from final 2000 library:

* `uaug_count > 0`: 451 / 2000.
* `uaug_count = 0`: 1549 / 2000.
* Top100 by heavy ensemble percentile/rank contained 1 uAUG-positive candidate.
* Top500 contained 8 uAUG-positive candidates.
* Upper ~50% heavy-ranked subset contained 26 uAUG-positive candidates among 927 heavy-ranked rows.
* Mean `heavy_ensemble_score`:
  * uAUG=0: 0.581325
  * uAUG=1: 0.316944
* Mean `robust_public_te_rank`:
  * uAUG=0: 0.6423
  * uAUG=1: 0.510599

Interpretation:

* uAUG is a strong negative feature in the current data/model behavior.
* High-ranked candidates are already strongly depleted for uAUG.
* The 451 uAUG-positive final rows are likely introduced mainly by quota/fill policy rather than by top-ranked model/evidence behavior.

uAUG by selection source:

```text
selection_source                   uAUG=0  uAUG=1
base_cand_diversity                   100       0
base_cand_exploratory                   4       8
base_cand_low_publicTE                 56      44
evidence_cand                        1006     232
fill_base_cand_clean_sequence          57      11
fill_evidence_cand                    326     156
```

PR3 planning implications:

* PR3 should start from `improved-v1.2`.
* PR3 should focus on selection policy rather than more PR2 bug-fixing.
* Current A-H quotas and evidence weights are heuristic and should be revisited.
* Blind fill should be replaced with explicit, auditable buckets.
* uAUG/uORF handling should become construct-aware:
  * exploitation/evidence buckets should strongly prefer or require uAUG-free candidates.
  * uAUG/uORF-positive candidates should be moved to designed controls or explicit exploration ladders rather than accidentally entering by fill logic.
* Restriction enzyme site checks should become soft/report-only for Gibson assembly unless they conflict with actual construct design.

Recommended PR3 branch:

```text
improved-v1.2-pr3-selection-policy
```

Recommended PR3 first step:

* Add reporting-only uAUG source audit and uAUG=0 dry-run selection before changing production output.

## 2026-06-04 - PR2 validation checkpoint: TSS extension and expressed-only TE labels

Branch:

* `improved-v1.1-pr2-tss-expression`

Latest commits:

* `d52210d` Add PR2 TSS extension and expressed-only TE labels
* `92b3279` Make PR2 selection expression gates group-specific

Status:

* Code pushed.
* `py_compile` passed by Codex for 02/03/10.
* Synthetic smoke test passed by Codex for group-specific 10 selection.
* Full raw-data validation was pending at this checkpoint; it later passed and was promoted to `improved-v1.2`.

Changed files:

```text
01_pipeline/scripts/02_tss_correction.py
01_pipeline/scripts/03_map_rna_ribo_public_te.py
01_pipeline/scripts/10_select_2000_cluster_diverse_library.py
docs/VALIDATION_PLAN.md
```

Summary:

* PR2 improves biological reliability of TSS correction and public TE labels.
* PR2 keeps PR1 08 -> 07 -> 10 heavy-score connection intact.
* PR2 keeps final selected CSV/FASTA filenames unchanged.
* PR2 changes final selection from global expression gate to group-specific expression gate.

Details:

1. `02_tss_correction.py`

* Added `trim_to_tss` / `extend_to_tss` behavior.
* Upstream TSS within `max_extend` can extend annotated 5'UTR.
* Added TSS correction mode and QC summary.

2. `03_map_rna_ribo_public_te.py`

* Added `is_expressed_public`.
* Added `expression_qc_reason`.
* TE/residual/`robust_public_te_rank` are calculated only for expressed rows.
* Non-expressed rows are treated as unreliable labels rather than forced low scores.

3. `10_select_2000_cluster_diverse_library.py`

* `base_cand` uses sequence QC only.
* `evidence_cand` uses `base_cand` + expression pass + `robust_public_te_rank`.
* A/B/C/D/E use `evidence_cand`.
* F/G/H can use clean `base_cand` where appropriate.
* J fill uses `evidence_cand` first, then `base_cand`.
* Added `selection_source` reporting.
* Added `evidence_candidate_pool_after_expression_TE_QC` reporting.

Validation required at checkpoint:

* Full pipeline rerun with real raw data.
* Check `selected_n` remains 2000.
* Check final library length remains 50-100.
* Check `max_per_seq_cluster` remains <= 2.
* Check `heavy_ensemble_score` present in final output.
* Check TSS correction summary.
* Check expressed-only label counts.
* Check `selection_source` counts.
* Rerun PR1.1 split disjointness check after PR2.

## 2026-06-04 - PR1.1 validated: split disjointness and final diversity reports

Branch:

* `improved-v1.1-pr1.1`
* merged into `improved-v1.1`

Commit:

* `b061ee7` Add PR1.1 split disjointness and final diversity reports

Changed files:

```text
01_pipeline/scripts/09_cluster_aware_classification_benchmark.py
01_pipeline/scripts/10_select_2000_cluster_diverse_library.py
```

Added outputs:

```text
06_modeling/tables/cluster_split_disjointness_check.csv
06_modeling/tables/final_library_gene_cluster_diversity_summary.txt
```

Validation:

* `gene_seq_cluster_split` passed for all four targets.
* Gene overlap = 0.
* Sequence cluster overlap = 0.
* `pass_required_for_split = True`.
* Final library `selected_n = 2000`.
* `n_unique_seq_clusters = 1943`.
* `max_per_seq_cluster = 2`.
* `n_unique_genes = 1850`.
* `max_per_gene = 4`.
* `heavy_ensemble_score` non-null = 2000 / 2000.
* length range = 50-100.

Interpretation:

* Final library is sequence-cluster-diverse.
* `seq_cluster_id` is a 5'UTR sequence-similarity cluster, not a gene cluster.
* Gene-level diversity is high but not hard-capped.

## 2026-06-03 - PR1 validated: cluster-aware heavy modeling connected to final selection

Branch:

* `improved-v1.1`

Commit:

* `4a29c08` Connect cluster-aware heavy modeling to final selection

Changed files:

```text
01_pipeline/scripts/run_00_full_final_pipeline.py
01_pipeline/scripts/07_heavy_rnafold_kmer6_automl.py
01_pipeline/scripts/10_select_2000_cluster_diverse_library.py
```

Changes:

* Run order changed so 08 Jaccard clustering runs before 07 heavy modeling.
* 07 now reads `tss_corrected_5utr_with_seq_clusters.csv` as primary input.
* 07 adds `seq_cluster_split` and `gene_seq_cluster_split`.
* 07 outputs `tss_corrected_5utr_with_seq_clusters_and_heavy_scores.csv`.
* 10 prioritizes heavy-score integrated input table.
* Final CSV/FASTA filenames remain unchanged.

Validation:

* `py_compile` passed for modified scripts.
* 07 created heavy-score integrated table.
* 09 benchmark completed.
* 10 selected final 2000-member library.
* `heavy_ensemble_score` non-null = 2000 / 2000.
* length range = 50-100.
* max per `seq_cluster_id` = 2.

Key benchmark result:

* `robust_public_te_rank` gene_seq_cluster_split RandomForest: AUC 0.671 / AP 0.582
* `protein_residual_rank` gene_seq_cluster_split ExtraTrees: AUC 0.696 / AP 0.630

Interpretation:

* Sequence features retain predictive signal even under strict gene + sequence-cluster split.
* `protein_residual_rank` is useful auxiliary evidence.
* `protein_abundance_rank` is weaker under strict split and should remain auxiliary.

## 2026-06-03 - Raw Data And LFS Handling

Summary:

* Raw public datasets were added using Git LFS where needed.
* Large NCBI mapping resources were reduced or filtered where possible.
* Code and raw data can be retrieved for reproducible local/company workstation runs.

Notes:

* Generated outputs should generally not be committed.
* Raw data may be kept via Git LFS or local/company storage depending on size/security constraints.
