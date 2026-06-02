# Data assets

`Heffner_minimal.tsv.gz` is a minimal gzip-compressed TSV derived from the public Heffner CHO proteomics supplementary table. It avoids NASCA/DRM issues caused by Office/CSV files in the company environment.

Use:

```bash
gunzip -c data_assets/Heffner_minimal.tsv.gz > 00_raw_data/05_cho_proteomics/Heffner_minimal.tsv
```
