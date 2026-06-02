from pathlib import Path
import argparse
import gzip
import json
import re
import numpy as np
import pandas as pd

BASE = Path.cwd()

PROT_DIR = BASE / "00_raw_data/05_cho_proteomics"
PROT_MINIMAL_CSV = PROT_DIR / "Heffner_minimal.csv"
PROT_MINIMAL_TSV = PROT_DIR / "Heffner_minimal.tsv"
PROT_XLSX = PROT_DIR / "Heffner_2020_CHO_hamster_proteomics_supp_table1.xlsx"

MAPDIR = PROT_DIR / "ncbi_gene_mapping"
GENE2ACC = MAPDIR / "gene2accession.gz"
GENEINFO = MAPDIR / "gene_info.gz"

OUT_CSV = PROT_DIR / "Heffner_2020_CHO_hamster_proteomics_ncbi_mapped_for_5utr.csv"
OUT_REPORT = PROT_DIR / "Heffner_2020_CHO_hamster_proteomics_ncbi_mapped_report.txt"
CFG_PATH = BASE / "01_pipeline/config/proteomics_config.json"

TAX_ID = 10029  # Cricetulus griseus


def norm_col(c):
    """Normalize column name for robust matching."""
    return re.sub(r"[^a-z0-9]+", "", str(c).lower())


def norm_symbol(x):
    if x is None or pd.isna(x):
        return None
    s = str(x).strip()
    if not s or s == "-":
        return None
    s = re.split(r"[;,|/\s]+", s)[0]
    s = re.sub(r"[^A-Za-z0-9_.-]", "", s)
    return s.upper() if s else None


def clean_accession(x):
    if x is None or pd.isna(x):
        return None
    s = str(x).strip()
    if not s or s.lower() in {"nan", "none", "-"}:
        return None
    # Remove Excel artifacts/spaces.
    s = s.replace('"', "").replace("'", "").strip()
    return s


def read_minimal_table():
    """Read Heffner_minimal.csv/tsv.

    Expected ideal header:
    Accession,GeneSymbol,PSMs,Peptides,UniquePeptides

    GeneSymbol/Peptides/UniquePeptides are optional, but Accession and PSMs are strongly recommended.
    """
    if PROT_MINIMAL_TSV.exists() and PROT_MINIMAL_TSV.stat().st_size > 0:
        path = PROT_MINIMAL_TSV
        forced_seps = ["\t", ","]
    elif PROT_MINIMAL_CSV.exists() and PROT_MINIMAL_CSV.stat().st_size > 0:
        path = PROT_MINIMAL_CSV
        forced_seps = [",", "\t"]
    else:
        raise SystemExit(
            "Minimal Heffner file not found.\n"
            f"Create one of:\n  {PROT_MINIMAL_CSV}\n  {PROT_MINIMAL_TSV}\n\n"
            "Required header example:\n"
            "Accession,GeneSymbol,PSMs,Peptides,UniquePeptides\n"
        )

    print("[3] FORCE load minimal Heffner file:", path)

    raw = path.read_bytes()
    print("  file size:", len(raw), "bytes")
    print("  first 80 bytes:", raw[:80])

    # If DRM file, fail clearly.
    if b"NASCA DRM FILE" in raw[:500] or b"DRM" in raw[:200]:
        raise SystemExit(
            "This minimal file still appears to be DRM-wrapped. Python cannot read NASCA DRM files.\n"
            "Create a brand-new plain text file outside DRM control, or ask IT/security for a DRM-free export."
        )

    last_err = None
    parsed = []
    for enc in ["utf-8-sig", "utf-8", "cp949", "euc-kr", "latin1"]:
        for sep in forced_seps + [None]:
            try:
                kwargs = dict(encoding=enc, engine="python", on_bad_lines="skip")
                if sep is not None:
                    kwargs["sep"] = sep
                else:
                    kwargs["sep"] = None
                df = pd.read_csv(path, **kwargs)
                df = df.dropna(axis=1, how="all").dropna(axis=0, how="all")
                df.columns = [str(c).strip().replace("\ufeff", "") for c in df.columns]
                parsed.append((df.shape[1], df.shape[0], enc, sep, df))
            except Exception as e:
                last_err = e

    if not parsed:
        raise RuntimeError(f"Could not read minimal file: {path}; last_error={last_err}")

    # Choose parse with the most columns, then rows.
    parsed.sort(key=lambda x: (x[0], x[1]), reverse=True)
    ncol, nrow, enc, sep, df = parsed[0]
    print(f"  parsed encoding={enc}, sep={repr(sep)}, shape={df.shape}")
    print("  columns:", list(df.columns))

    # If first row accidentally contains header, promote it.
    normalized = {norm_col(c): c for c in df.columns}
    if not any(k in normalized for k in ["accession", "proteinaccession", "acc"]):
        # Try to find a row with Accession/PSM header and reread with skiprows.
        text = raw.decode(enc, errors="replace")
        lines = text.splitlines()
        header_idx = None
        for i, line in enumerate(lines[:300]):
            low = line.lower()
            if "accession" in low and ("psm" in low or "peptide" in low or "gene" in low):
                header_idx = i
                break
        if header_idx is not None:
            sep2 = "\t" if lines[header_idx].count("\t") > lines[header_idx].count(",") else ","
            df = pd.read_csv(path, encoding=enc, sep=sep2, engine="python", skiprows=header_idx, on_bad_lines="skip")
            df = df.dropna(axis=1, how="all").dropna(axis=0, how="all")
            df.columns = [str(c).strip().replace("\ufeff", "") for c in df.columns]
            print(f"  reparsed with header line {header_idx+1}, sep={repr(sep2)}, shape={df.shape}")
            print("  columns:", list(df.columns))

    return df


def find_col(df, candidates, required=False):
    norm_map = {norm_col(c): c for c in df.columns}

    # direct normalized candidates
    for cand in candidates:
        key = norm_col(cand)
        if key in norm_map:
            return norm_map[key]

    # contains candidates
    for c in df.columns:
        nc = norm_col(c)
        for cand in candidates:
            key = norm_col(cand)
            if key and key in nc:
                return c

    if required:
        raise KeyError(
            f"Could not find required column among candidates={candidates}\n"
            f"Available columns={list(df.columns)}\n"
            "Make sure Heffner_minimal.csv first row is exactly:\n"
            "Accession,GeneSymbol,PSMs,Peptides,UniquePeptides"
        )
    return None


def read_gene_info():
    if not GENEINFO.exists():
        print("[WARN] gene_info.gz missing; will rely on GeneSymbol if present.")
        return pd.DataFrame(columns=["gene_id_key", "gene_symbol_key", "Symbol", "description"])

    cols = [
        "tax_id", "GeneID", "Symbol", "LocusTag", "Synonyms", "dbXrefs", "chromosome",
        "map_location", "description", "type_of_gene", "Symbol_from_nomenclature_authority",
        "Full_name_from_nomenclature_authority", "Nomenclature_status", "Other_designations",
        "Modification_date", "Feature_type"
    ]
    use = ["tax_id", "GeneID", "Symbol", "Synonyms", "description"]
    print("[1] Load gene_info.gz for tax_id", TAX_ID)
    chunks = []
    for chunk in pd.read_csv(GENEINFO, sep="\t", names=cols, comment="#", dtype=str, chunksize=500000):
        sub = chunk[chunk["tax_id"] == str(TAX_ID)][use].copy()
        if len(sub):
            chunks.append(sub)
    if not chunks:
        print("[WARN] No Cricetulus griseus rows found in gene_info.gz")
        return pd.DataFrame(columns=["gene_id_key", "gene_symbol_key", "Symbol", "description"])

    gi = pd.concat(chunks, ignore_index=True)
    gi["gene_id_key"] = gi["GeneID"].astype(str)
    gi["gene_symbol_key"] = gi["Symbol"].apply(norm_symbol)
    print("  gene_info rows:", len(gi))
    return gi


def read_gene2accession():
    if not GENE2ACC.exists():
        print("[WARN] gene2accession.gz missing; will rely on GeneSymbol if present.")
        return pd.DataFrame(columns=[
            "protein_gi", "protein_accession.version", "gene_id_key", "gene_symbol_key_gene2acc"
        ])

    cols = [
        "tax_id", "GeneID", "status", "RNA_nucleotide_accession.version",
        "RNA_nucleotide_gi", "protein_accession.version", "protein_gi",
        "genomic_nucleotide_accession.version", "genomic_nucleotide_gi",
        "start_position_on_the_genomic_accession", "end_position_on_the_genomic_accession",
        "orientation", "assembly", "mature_peptide_accession.version",
        "mature_peptide_gi", "Symbol"
    ]
    use = ["tax_id", "GeneID", "protein_accession.version", "protein_gi", "Symbol"]
    print("[2] Load gene2accession.gz for tax_id", TAX_ID)
    chunks = []
    for chunk in pd.read_csv(GENE2ACC, sep="\t", names=cols, comment="#", dtype=str, chunksize=1000000):
        sub = chunk[chunk["tax_id"] == str(TAX_ID)][use].copy()
        if len(sub):
            chunks.append(sub)
    if not chunks:
        print("[WARN] No Cricetulus griseus rows found in gene2accession.gz")
        return pd.DataFrame(columns=[
            "protein_gi", "protein_accession.version", "gene_id_key", "gene_symbol_key_gene2acc"
        ])

    g2a = pd.concat(chunks, ignore_index=True)
    g2a["gene_id_key"] = g2a["GeneID"].astype(str)
    g2a["gene_symbol_key_gene2acc"] = g2a["Symbol"].apply(norm_symbol)
    print("  gene2accession rows:", len(g2a))
    return g2a


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--show-unmapped", action="store_true")
    args = parser.parse_args()

    gi = read_gene_info()
    g2a = read_gene2accession()
    raw = read_minimal_table()

    acc_col = find_col(raw, ["Accession", "ProteinAccession", "Protein accession", "Acc"], required=True)
    gene_col = find_col(raw, ["GeneSymbol", "Gene Symbol", "Symbol", "Gene", "gene_symbol"], required=False)
    psm_col = find_col(raw, ["PSMs", "PSM", "PSMCount", "PSM Count", "TotalPSMs", "Total PSMs"], required=True)
    pep_col = find_col(raw, ["Peptides", "Peptide", "PeptideCount", "Peptide Count", "TotalPeptides"], required=False)
    uniq_col = find_col(raw, ["UniquePeptides", "Unique Peptides", "UniquePeptide", "Unique peptide count"], required=False)
    desc_col = find_col(raw, ["Description", "ProteinDescription", "Protein Description"], required=False)

    prot = pd.DataFrame()
    prot["Accession"] = raw[acc_col].map(clean_accession)
    prot["GeneSymbol"] = raw[gene_col].map(norm_symbol) if gene_col else None
    prot["Description"] = raw[desc_col].astype(str) if desc_col else prot["GeneSymbol"].fillna("").astype(str)
    prot["Σ# PSMs"] = pd.to_numeric(raw[psm_col], errors="coerce")
    prot["Σ# Peptides"] = pd.to_numeric(raw[pep_col], errors="coerce") if pep_col else prot["Σ# PSMs"]
    prot["Σ# Unique Peptides"] = pd.to_numeric(raw[uniq_col], errors="coerce") if uniq_col else prot["Σ# Peptides"]
    prot["ΣCoverage"] = 0
    prot["source_sheet"] = "minimal_manual"

    prot = prot[prot["Accession"].notna() | prot["GeneSymbol"].notna()].copy()
    prot = prot[prot["Σ# PSMs"].notna()].copy()

    print("[4] Minimal protein rows usable:", len(prot))
    print(prot.head(5).to_string(index=False))

    # accession mapping
    prot["protein_gi_key"] = pd.to_numeric(prot["Accession"], errors="coerce").astype("Int64").astype(str)
    prot.loc[prot["protein_gi_key"] == "<NA>", "protein_gi_key"] = None
    prot["protein_accession_key"] = prot["Accession"].astype(str)
    prot["protein_accession_noversion"] = prot["protein_accession_key"].str.replace(r"\.\d+$", "", regex=True)

    mapped = prot.copy()
    mapped["gene_id_key"] = None
    mapped["gene_symbol_key_gene2acc"] = None

    if len(g2a):
        # map numeric protein_gi
        g2a_gi = g2a[["protein_gi", "gene_id_key", "gene_symbol_key_gene2acc"]].dropna().drop_duplicates("protein_gi")
        m1 = mapped.merge(g2a_gi, left_on="protein_gi_key", right_on="protein_gi", how="left", suffixes=("", "_m"))
        mapped["gene_id_key"] = m1["gene_id_key_m"].combine_first(mapped["gene_id_key"])
        mapped["gene_symbol_key_gene2acc"] = m1["gene_symbol_key_gene2acc_m"].combine_first(mapped["gene_symbol_key_gene2acc"])

        # map protein accession.version
        need = mapped["gene_id_key"].isna()
        if need.any():
            g2a_acc = g2a[["protein_accession.version", "gene_id_key", "gene_symbol_key_gene2acc"]].dropna().drop_duplicates("protein_accession.version")
            fb = mapped.loc[need].merge(g2a_acc, left_on="protein_accession_key", right_on="protein_accession.version", how="left")
            mapped.loc[need, "gene_id_key"] = fb["gene_id_key_y"].values
            mapped.loc[need, "gene_symbol_key_gene2acc"] = fb["gene_symbol_key_gene2acc_y"].values

        # map protein accession no version
        need = mapped["gene_id_key"].isna()
        if need.any():
            g2a_acc2 = g2a[["protein_accession.version", "gene_id_key", "gene_symbol_key_gene2acc"]].dropna().copy()
            g2a_acc2["protein_accession_noversion"] = g2a_acc2["protein_accession.version"].astype(str).str.replace(r"\.\d+$", "", regex=True)
            g2a_acc2 = g2a_acc2.drop_duplicates("protein_accession_noversion")
            fb = mapped.loc[need].merge(g2a_acc2, on="protein_accession_noversion", how="left")
            mapped.loc[need, "gene_id_key"] = fb["gene_id_key_y"].values
            mapped.loc[need, "gene_symbol_key_gene2acc"] = fb["gene_symbol_key_gene2acc_y"].values

    if len(gi) and "gene_id_key" in gi.columns:
        mapped = mapped.merge(gi[["gene_id_key", "gene_symbol_key", "Symbol", "description"]], on="gene_id_key", how="left")
    else:
        mapped["gene_symbol_key"] = None
        mapped["Symbol"] = None
        mapped["description"] = None

    # fallback to manual GeneSymbol if gene mapping failed
    mapped["gene_symbol_key"] = mapped["gene_symbol_key"].fillna(mapped["gene_symbol_key_gene2acc"]).fillna(mapped["GeneSymbol"])
    mapped["gene_symbol"] = mapped["Symbol"].fillna(mapped["gene_symbol_key"])

    # If no gene_id but symbol exists, create pseudo gene_id. Downstream can still map by gene_symbol.
    no_id_has_symbol = mapped["gene_id_key"].isna() & mapped["gene_symbol_key"].notna()
    mapped.loc[no_id_has_symbol, "gene_id_key"] = "SYMBOL_" + mapped.loc[no_id_has_symbol, "gene_symbol_key"].astype(str)

    proxy_cols = ["Σ# PSMs", "Σ# Peptides", "Σ# Unique Peptides"]
    tmp = mapped[proxy_cols].apply(pd.to_numeric, errors="coerce").fillna(0)
    mapped["abundance_proxy"] = np.log1p(tmp.mean(axis=1))

    out = mapped[mapped["gene_symbol_key"].notna() & mapped["gene_id_key"].notna()].copy()

    if len(out) == 0:
        OUT_REPORT.write_text(
            "0 mapped rows from minimal Heffner file.\n"
            f"input rows: {len(prot)}\n"
            f"columns raw: {list(raw.columns)}\n"
            f"selected acc_col={acc_col}, gene_col={gene_col}, psm_col={psm_col}, pep_col={pep_col}, uniq_col={uniq_col}\n"
            "\nFirst rows:\n"
            + prot.head(30).to_string(index=False),
            encoding="utf-8"
        )
        raise SystemExit(f"0 usable mapped rows. See {OUT_REPORT}")

    keep = [
        "gene_id_key", "gene_symbol_key", "gene_symbol", "abundance_proxy",
        "Description", "Accession", "source_sheet",
        "Σ# PSMs", "Σ# Peptides", "Σ# Unique Peptides", "ΣCoverage", "description"
    ]

    agg = {
        "gene_symbol_key": "first",
        "gene_symbol": "first",
        "abundance_proxy": "mean",
        "Description": "first",
        "Accession": "first",
        "source_sheet": lambda x: ";".join(sorted(set(map(str, x)))),
        "Σ# PSMs": "mean",
        "Σ# Peptides": "mean",
        "Σ# Unique Peptides": "mean",
        "ΣCoverage": "mean",
        "description": "first",
    }

    collapsed = out[keep].groupby("gene_id_key", as_index=False).agg(agg)
    OUT_CSV.parent.mkdir(parents=True, exist_ok=True)
    collapsed.to_csv(OUT_CSV, index=False, encoding="utf-8")

    # update config
    if CFG_PATH.exists():
        try:
            cfg = json.loads(CFG_PATH.read_text(encoding="utf-8"))
        except Exception:
            cfg = {}
    else:
        cfg = {}

    cfg.update({
        "proteomics_path": "00_raw_data/05_cho_proteomics/Heffner_2020_CHO_hamster_proteomics_ncbi_mapped_for_5utr.csv",
        "gene_symbol_column": "gene_symbol_key",
        "gene_id_column": "gene_id_key",
        "abundance_columns": ["abundance_proxy"],
        "protein_abundance_transform": "none",
        "protein_aggregation": "mean",
    })
    cfg.setdefault("sheet_name", 0)
    cfg.setdefault("multi_omics_score_weights", {
        "robust_public_te_rank": 0.50,
        "protein_residual_rank": 0.20,
        "protein_abundance_rank": 0.15,
        "day_consensus_TE_rank": 0.10,
        "tss_confidence_score": 0.05
    })
    cfg.setdefault("candidate_filter", {
        "training_min_length": 40,
        "training_max_length": 200,
        "selection_min_length": 50,
        "selection_max_length": 100,
        "gc_min": 0.30,
        "gc_max": 0.75,
        "uaug_max": 1
    })
    CFG_PATH.parent.mkdir(parents=True, exist_ok=True)
    CFG_PATH.write_text(json.dumps(cfg, indent=2), encoding="utf-8")

    report = [
        "FORCED minimal Heffner proteomics mapping report",
        "=" * 90,
        f"input_file: {PROT_MINIMAL_TSV if PROT_MINIMAL_TSV.exists() else PROT_MINIMAL_CSV}",
        f"raw_rows: {len(raw)}",
        f"usable_protein_rows: {len(prot)}",
        f"rows_with_gene_symbol_or_mapping: {len(out)}",
        f"unique_gene_ids_or_pseudo_ids: {collapsed['gene_id_key'].nunique()}",
        f"unique_gene_symbols: {collapsed['gene_symbol_key'].nunique()}",
        f"selected columns: acc={acc_col}, gene={gene_col}, psm={psm_col}, peptides={pep_col}, unique={uniq_col}, desc={desc_col}",
        f"output_csv: {OUT_CSV}",
        f"updated_config: {CFG_PATH}",
        "",
        "[Top 30 by abundance_proxy]",
        collapsed.sort_values("abundance_proxy", ascending=False).head(30)[
            ["gene_id_key", "gene_symbol_key", "gene_symbol", "abundance_proxy", "Description", "Accession"]
        ].to_string(index=False),
    ]
    OUT_REPORT.write_text("\n".join(report), encoding="utf-8")

    print("[SAVED]", OUT_CSV)
    print("[SAVED]", OUT_REPORT)
    print("[UPDATED]", CFG_PATH)
    print("\n".join(report[:12]))

    if args.show_unmapped:
        unmapped = mapped[mapped["gene_symbol_key"].isna()][["Accession", "GeneSymbol", "Description"]].head(100)
        print("\n[Unusable examples]")
        print(unmapped.to_string(index=False))


if __name__ == "__main__":
    main()
