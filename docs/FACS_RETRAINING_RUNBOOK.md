# FACS-guided retraining runbook

## Purpose

This follow-up stage connects the final v1.4 selected 5'UTR table to the experimental
SortSeq phenotype table by the UTR sequence itself. It answers three separate questions:

1. How the measured phenotypes and high-expression hits are distributed across the exact
   `library_group` categories used in the 2,000-UTR library.
2. How well each original public-TE/model/multiomics score agrees with the new FACS result.
3. How well a FACS-trained model predicts held-out genes and sequence-similarity clusters.

It does not alter or rerun steps 00-11. Run it only after both the v1.4 selected-library CSV
and SortSeq `utr_results_full.tsv` are present on the offline Linux server.

## Input contract

### Selected library

Default path:

```text
07_library_design/tables/selected_2000_50_100bp_cluster_diverse_evidence_balanced_library.csv
```

The join key is normalized `utr5_sequence_tss_corrected`. The script preserves the actual
`library_group` strings found in this file; it does not replace them with expected quota names.
This is important because the v1.4 K-refill groups can make the realized eight-category counts
different from the primary A-H quota declaration in the selection script.

### SortSeq phenotype table

Use the full table from the separate `5UTR_SortSeq` pipeline:

```text
/data/user/MCET03/03_NGS/02_5UTR_sorting/results/sortseq/utr_results_full.tsv
```

SortSeq defines `variant_id` from the LibraryQC `sequence_key`, and also normally retains a
`sequence` column. The retraining script auto-detects either and performs a normalized
sequence-to-sequence one-to-one join. A match rate below 95% stops the run.

## Experimental endpoints

Both endpoints are retained rather than forcing them into one score:

| Endpoint | Role | Meaning |
|---|---|---|
| `high15_probability` | Primary | Estimated probability of bin1 or bin2 within the mCherry+/GFP- target gate |
| `expected_bin_score` | Secondary | Whole-distribution weighted average using bin1-bin6 weights 6,5,4,3,2,1 |

Default eligibility is `top15_read_support_pass=TRUE` together with total six-bin count at least
200. This carries forward the SortSeq High15 read-support rule (including high-bin support). The
unsorted sample is not used as the denominator of either fluorescence endpoint.

## Run on the offline server

From the root of the `5UTR_Engineering` branch:

```bash
FACS_N_JOBS=64 bash RUN_FACS_RETRAINING.sh \
  /data/user/MCET03/03_NGS/02_5UTR_sorting/results/sortseq/utr_results_full.tsv
```

The wrapper uses the v1.4 selected-library file in the current repository. To override any
setting, append normal Python options, for example:

```bash
FACS_N_JOBS=64 bash RUN_FACS_RETRAINING.sh \
  /data/user/MCET03/03_NGS/02_5UTR_sorting/results/sortseq/utr_results_full.tsv \
  --n-estimators 600 --min-match-rate 0.99
```

For a fast input/merge audit without model fitting:

```bash
bash RUN_FACS_RETRAINING.sh \
  /data/user/MCET03/03_NGS/02_5UTR_sorting/results/sortseq/utr_results_full.tsv \
  --skip-modeling --no-plots
```

## Main outputs

All outputs are written to `08_reports/facs_retraining/`.

| File | Question answered |
|---|---|
| `FACS_RETRAINING_SUMMARY.md` | One-page run summary and interpretation boundary |
| `run_manifest.json` | Input hashes, join columns, match rate, group counts, parameters and model paths |
| `tables/category_distribution_and_scores.csv` | Exact eight-group sizes, eligible counts, High15 and expected-score summaries |
| `tables/category_by_facs_high15_tier.csv` | Within-category distribution across Top 5%, 5-15%, 15-25%, 25-50% and bottom 50% |
| `tables/category_top_hit_enrichment.csv` | Group-specific Top 5/10/15/20% hit rates and observed/expected enrichment |
| `tables/original_score_vs_facs_correlations.csv` | Spearman, Pearson, Kendall and percentile gap for each old score/end point pair |
| `tables/original_score_vs_facs_topn_overlap.csv` | Top 20/50/100/200 overlap and random-overlap expectation |
| `tables/facs_model_oof_metrics.csv` | Leakage-aware OOF performance for sequence-only, prior-only and hybrid models |
| `tables/facs_model_oof_predictions.csv` | Variant-level observed values, fold, leakage group and every OOF prediction |
| `models/*.joblib` | Final models fit to all eligible variants after OOF evaluation |

The category enrichment table is preferable to raw hit counts because the eight library groups
have different sizes. A group with many Top15 hits can still be depleted after dividing by its
number of eligible UTRs.

## Model validation design

The script makes connected components using both `seq_cluster_id` and gene identity. If two UTRs
share either a gene or a sequence cluster, they are assigned to the same cross-validation fold.
This prevents a near-duplicate or same-gene candidate from appearing in training while its close
counterpart is treated as held-out test data.

Three ExtraTrees regressors are compared for both FACS endpoints:

- `sequence_only`: length/composition/uORF/homopolymer and 2-6-mer features.
- `prior_only`: original prediction/evidence scores present in the selected-library table.
- `hybrid`: sequence features plus the original scores.

Only OOF predictions are used for the reported Spearman, MAE, R2, Top15 AUROC/AP and Top-N
overlap. The saved final models are refit on all eligible rows after this evaluation and therefore
must not be described using their in-sample performance.

## Interpretation boundary

This dataset comes from one pooled SortSeq experiment. Variant-level read bootstrap and OOF model
validation do not replace biological replicate uncertainty. Results support statements about
relative steady-state mCherry fluorescence within this selected 2,000-UTR design space. They do
not by themselves establish absolute translation rate, a population-wide causal effect, or
generalization to arbitrary CHO 5'UTRs.
