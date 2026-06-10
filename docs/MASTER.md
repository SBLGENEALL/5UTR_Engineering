# 5UTR Engineering MASTER

Last updated: 2026-06-10

## 1. Current Project State

Repository:

* `SBLGENEALL/5UTR_Engineering`

Stable branch:

* `main`
* Initial stable numbered pipeline release.
* Do not replace `main` until PR2/PR3 and release documentation are complete.

Current validated release branch:

* `v1.3`
* Defined from validated PR3-2 uAUG0 production selection.
* PR3-2 workstation validation passed.
* This is the clean baseline to preserve before starting v1.4 development.

Validated v1.3 final library QC:

* `selected_n = 2000`
* `requested_n = 2000`
* `uaug_positive_n = 0`
* `uaug0_policy_pass = True`
* `n_unique_seq_clusters = 1937`
* `max_per_seq_cluster = 2`
* `cluster_cap = 2`
* `cluster_cap_pass = True`
* `gene_key = gene_name`
* `n_unique_genes = 1896`
* `max_per_gene = 3`
* `gene_cap = 4`
* `gene_cap_pass = True`
* `mean_heavy_ensemble_score = 0.5811515`
* `mean_robust_public_te_rank = 0.642380`

Release handling:

* Do not merge `v1.3` into `main` as part of this QC staging commit.
* Do not include raw data, caches, logs, or intermediate outputs.
* Preserve only the validated final-library QC summary in Git.

Previous validated baseline:

* `improved-v1.1`
* PR1 + PR1.1 validated.
* The 08 Jaccard sequence clustering -> 07 heavy RNAfold/k-mer modeling -> 09 cluster-aware benchmark -> 10 final library selection connection was validated.

Validated PR2 branch:

* `improved-v1.1-pr2-tss-expression`
* Purpose: TSS correction + expressed-only public TE labels + group-specific final selection expression gate.
* Status: full raw-data validation passed; suitable to promote as `improved-v1.2` working baseline.

Current operating rule:

* Treat this `docs/MASTER.md` and `docs/CHANGELOG.md` as the source of truth for this branch family.
* Historical chat context is secondary.
* Workspace chats are disposable; validated decisions should be written back to GitHub docs.

## 2. Current Validated Baseline: improved-v1.2

PR2 goals:

* Improve TSS correction using atlas-supported trim/extend behavior.
* Improve public TE label reliability by computing TE/residual labels only for expressed rows.
* Preserve PR1 08 -> 07 -> 10 heavy-score connection.
* Keep final selected CSV/FASTA filenames unchanged.
* Apply expression gating to evidence-supported groups while allowing clean exploratory/diversity/negative-control candidates where appropriate.

Validated PR2 results from company workstation:

* Final selected library generated successfully.
* `selected_n = 2000`.
* CSV line count = 2001 including header.
* Length mean = 73.527.
* Length min/max = 50/100.
* `max_per_cluster_primary = 1`.
* `heavy_ensemble_score` and `heavy_ensemble_rank` are present in the final selected library.

TSS correction validation:

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

Interpretation:

* PR2 improved the strict TE classification signal relative to PR1.
* The likely driver is reduced sequence/label noise from atlas-based TSS correction and expressed-only TE labels.
* ExtraTrees is the best observed classifier in this checkpoint, but final selection uses the integrated heavy/model/evidence columns rather than treating a single classifier as the only criterion.

Leakage/disjointness validation:

* For the strict `gene_seq_cluster_split`, overlap/leakage was 0 for all four targets checked:
  * TE
  * multi-omics
  * protein abundance
  * protein residual
* Non-strict split modes can show overlaps by design; the required strict split passed.

Final selection source counts:

```text
evidence_cand                 1238
fill_base_cand_clean_sequence  482
base_cand_diversity            100
base_cand_low_publicTE         100
fill_evidence_cand              68
base_cand_exploratory           12
```

Interpretation:

* The final 2,000-member PR2 library is evidence-driven, not simply the top 2,000 by model score.
* Most candidates come from expressed public TE evidence candidates.
* A substantial clean sequence fill component exists because the primary selection enforces cluster diversity.

## 3. Pipeline Order

Expected full pipeline order:

```text
00_check_inputs.py
01_build_utr_database.py
02_tss_correction.py
03_map_rna_ribo_public_te.py
04_preprocess_heffner_proteomics.py
05_integrate_proteomics_multiomics.py
06_plot_multiomics_distributions.py
08_jaccard_sequence_cluster_qc.py
07_heavy_rnafold_kmer6_automl.py
09_cluster_aware_classification_benchmark.py
10_select_2000_cluster_diverse_library.py
```

Important:

* 08 must run before 07 because 07 uses `seq_cluster_id`.
* 07 generates `heavy_ensemble_score` and the integrated heavy-score table.
* 10 should preferentially read:

```text
04_te_labeling/tables/tss_corrected_5utr_with_seq_clusters_and_heavy_scores.csv
```

## 4. PR2 Selection Logic Snapshot

The current PR2 selector first creates `base_cand` with sequence QC:

* length 50-100
* GC 0.30-0.75
* `uaug_count <= 1`
* no forbidden restriction-enzyme sites under the current PR2 logic
* non-empty sequence
* no `N`
* duplicate sequence removal

Then it creates `evidence_cand`:

* `base_cand`
* `is_expressed_public = True`
* non-null `robust_public_te_rank`

The current evidence score weights are heuristic, not yet performance-optimized:

```text
robust_public_te_rank      42%
day_consensus_TE_rank      16%
protein_residual_rank      12%
protein_abundance_rank     10%
multi_omics_utr_rank       10%
model_score                 7%
tss_confidence_score        3%
```

The original A-H quotas are also heuristic:

```text
A_publicTE_high_confidence          500
B_TE_model_classifier_supported     300
C_protein_abundance_supported       250
D_protein_residual_supported        250
E_multiomics_consensus_high         250
F_sequence_diverse_exploratory      200
G_length_GC_uAUG_diversity          100
H_low_signal_negative_controls      100
```

Because `max_per_cluster_primary = 1`, actual filled counts can differ strongly from the requested quotas.

## 5. PR3 Direction

PR3 should start from `improved-v1.2`.

Recommended branch name:

```text
improved-v1.2-pr3-selection-policy
```

PR3 goals:

* Make uAUG/uORF screening stricter for final CHO vector candidates.
* Treat restriction enzyme sites as soft warning or report-only because the planned library construction uses Gibson assembly, unless a site conflicts with the actual construct design.
* Reconsider the A-H fixed quota policy.
* Consider model-performance-informed weighting rather than only biology-intuition weights.
* Consider outputting alternative final libraries for comparison:
  * evidence-balanced
  * model-prioritized
  * conservative construct-ready

PR3 notes:

* PR2 still allows `uaug_count <= 1`; it does not fully remove uAUG/uORF risk.
* PR3 should distinguish:
  * hard fail: strong uORF/uAUG risk, invalid sequence, serious construct incompatibility
  * soft warning: enzyme sites under Gibson-compatible design, moderate motif concerns
  * report-only: descriptive features such as length/GC/structure bins

## 6. Git Strategy

`main`:

* Stable release only.

`improved-v1.1`:

* PR1 + PR1.1 validated historical baseline.

`improved-v1.2`:

* PR2 validated working baseline.
* Recommended base for PR3.

`improved-v1.2-pr3-selection-policy`:

* Planned PR3 branch.

Merge policy:

* Do not merge `improved-v1.2` into `main` until PR3 selection policy and release documentation are complete.
* Generated outputs should not be committed unless explicitly chosen as small release artifacts.
