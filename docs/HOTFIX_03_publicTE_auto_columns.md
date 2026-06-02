# HOTFIX: 03_map_rna_ribo_robust_public_te.py

## Issue

The final ZIP release version of `01_pipeline/scripts/03_map_rna_ribo_robust_public_te.py` can fail when `project_config.json` contains:

```json
"rna_day3": "auto",
"rna_day6": "auto",
"ribo_day3": "auto",
"ribo_day6": "auto"
```

The error is:

```text
TypeError: can only concatenate list (not "str") to list
```

## Cause

The script concatenates sample groups as lists, but the config value `auto` is read as a string.

## Fix

Replace `01_pipeline/scripts/03_map_rna_ribo_robust_public_te.py` with the patched version that resolves `auto` to:

```text
RNA day3:  s01, s02, s03
RNA day6:  s07, s08, s09
Ribo day3: s04, s05, s06
Ribo day6: s10, s11, s12
```

The patched script also prints the resolved columns before calculation.

## Expected successful log

```text
[INFO] RNA gene column: ...
[INFO] Ribo gene column: ...
[INFO] RNA day3 columns: ['s01', 's02', 's03']
[INFO] RNA day6 columns: ['s07', 's08', 's09']
[INFO] Ribo day3 columns: ['s04', 's05', 's06']
[INFO] Ribo day6 columns: ['s10', 's11', 's12']
[SAVED] 04_te_labeling/tables/tss_corrected_5utr_robust_public_te_labels.csv
[SAVED] 04_te_labeling/tables/tss_corrected_5utr_50_100bp_training_ready.csv
```

## Temporary workaround

Instead of using `auto`, edit `01_pipeline/config/project_config.json`:

```json
"rna_day3": ["s01", "s02", "s03"],
"rna_day6": ["s07", "s08", "s09"],
"ribo_day3": ["s04", "s05", "s06"],
"ribo_day6": ["s10", "s11", "s12"]
```

The script-level fix is preferred for the final team release.
