# GitHub release guide

## Recommended GitHub contents

Commit the clean numbered pipeline, docs, config templates, and `data_assets/Heffner_minimal.tsv.gz`.
Do not commit large raw data or generated intermediate results.

```bash
git add README.md RUN_FINAL_MAIN.sh RUN_FINAL_MAIN.bat .gitignore   docs data_assets 01_pipeline/config 01_pipeline/scripts tools   00_raw_data 02_utr_database 03_tss_correction 04_te_labeling   05_feature_extraction 06_modeling 07_library_design 08_reports

git commit -m "Finalize numbered CHO 5UTR pipeline with required TSS atlas"
git tag -a v1.0.0 -m "CHO 5UTR numbered final pipeline"
git push origin main --tags
```

## Raw data policy

Raw data files are not committed. Only folder placeholders are included.
Team users place raw data according to `docs/INPUTS_REQUIRED_TSS_RAWDATA.md` and run:

```bash
bash RUN_FINAL_MAIN.sh
```
