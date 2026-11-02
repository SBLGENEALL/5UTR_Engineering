# 08_TOP_candidate_analysis

Separate exploratory analysis; **no FACS retraining**. Based on v1.4-pr2.

## Inputs
- `hits.txt`: one UTR ID per line, e.g. UTR0552
- `--library`: final 2,000 library CSV/TSV (expected original output `07_library_design/tables/selected_2000_50_100bp_cluster_diverse_evidence_balanced_library.csv`)
- Optional `--validation`: CSV/TSV with one row per tested UTR, columns `UTR_ID`, `batch`, `relative_gMFI`. Include all 60 tested candidates for primary hit-vs-non-hit analysis.

## Run from repository root
```bash
Rscript 08_TOP_candidate_analysis/analyze_top_candidates.R --library /absolute/path/to/selected_2000.csv --hits 08_TOP_candidate_analysis/hits.txt --out 08_TOP_candidate_analysis/results
Rscript 08_TOP_candidate_analysis/analyze_top_candidates.R --library /absolute/path/to/selected_2000.csv --hits 08_TOP_candidate_analysis/hits.txt --validation /absolute/path/to/validation.csv --out 08_TOP_candidate_analysis/results
```

The original selection script uses `library_index` and `library_group`, but **UTR0552 = library_index 552 is not proven**. If the final table has an explicit `utr_id` or `UTR_ID`, it is used automatically. Otherwise inspect the ID mapping and use `--id-mode index` **only after confirming that correspondence**. The program intentionally stops rather than silently misassigning sequences.

GC and length can be computed from sequence if absent. MFE is used only if an actual MFE column exists; it is not invented. The original repository's generated results are generally local and may not be tracked in GitHub.

## Outputs
CSV: hit_feature_table, feature_tests, category_enrichment, expression_correlations, validation_checked (when supplied), hit_ids_used.
PDF: feature_boxplots, category_composition. TXT: analysis_summary.

Statistics: Wilcoxon rank sum for numeric features, Fisher exact for category, Spearman for expression correlations, Holm correction within comparison. This is **exploratory**, not evidence of causal 5'UTR determinants. Batch-adjusted inference is not implemented. Hit selection is relative top-7 within each batch, not an independent random sample.

R dependency: base R only.
