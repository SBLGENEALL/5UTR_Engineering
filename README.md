# 5′UTR Engineering Pipeline

This repository is the final team-distribution repository for the 5′UTR engineering project.

## Purpose

This pipeline organizes 5′UTR candidate ranking into a reproducible workflow:

1. Validate input tables
2. Build a training table
3. Remove highly similar sequences by k-mer/Jaccard similarity
4. Train a TE percentile-rank prediction model
5. Score candidate 5′UTR sequences
6. Select top candidates for synthesis/validation
7. Generate summary figures

## Repository structure

```text
.
├── config/                 # YAML configuration
├── data/
│   ├── input/              # user-provided training table templates
│   └── candidates/         # candidate sequence templates
├── docs/                   # pipeline documentation
├── scripts/                # executable analysis scripts
├── results/                # generated outputs; ignored by git except .gitkeep
├── archive_do_not_run/     # old or deprecated material
├── run_pipeline.py         # one-command runner
├── QUICK_START.md
├── VERSION
└── CHANGELOG.md
```

## Quick start

```bash
python scripts/00_check_environment.py
python run_pipeline.py --config config/config.yaml
```

## Main output

The most important output is:

```text
results/06_selection/top_candidates.csv
```

This file contains the final ranked 5′UTR candidates for experimental follow-up.

## Version

Current release: `v1.0.0`
