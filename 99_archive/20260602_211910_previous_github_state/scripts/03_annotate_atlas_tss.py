import argparse
import gzip
from pathlib import Path

import pandas as pd


def open_text(path: Path):
    return gzip.open(path, "rt", encoding="utf-8", errors="ignore") if path.suffix == ".gz" else path.open("r", encoding="utf-8", errors="ignore")


def read_tss_bed(path: Path) -> pd.DataFrame:
    cols = ["chromosome", "bed_start", "bed_end", "tss_id", "tss_score", "strand"]
    df = pd.read_csv(path, sep="\t", names=cols, compression="infer")
    df["tss_pos_1based"] = df["bed_start"].astype(int) + 1
    parts = df["tss_id"].astype(str).str.extract(r"^(?P<promoter_id>[^@]+)@(?P<gene_symbol>.+)_(?P<transcript_id>[A-Z]{2}_\d+\.\d+)$")
    df = pd.concat([df, parts], axis=1)
    return df


def read_tss_meta(path: Path) -> pd.DataFrame:
    meta = pd.read_csv(path, sep="\t", compression="infer")
    meta = meta.rename(columns={c: c.strip().replace(" ", "_") for c in meta.columns})
    first_col = meta.columns[0]
    if first_col != "tss_id":
        meta = meta.rename(columns={first_col: "tss_id"})
    return meta


def pick_best_tss(df: pd.DataFrame) -> pd.DataFrame:
    # Keep the strongest atlas TSS per transcript.
    # Score represents atlas support intensity; higher is preferred.
    df = df.copy()
    df["tss_score"] = pd.to_numeric(df["tss_score"], errors="coerce").fillna(0)
    df = df.sort_values(["transcript_id", "tss_score"], ascending=[True, False])
    return df.drop_duplicates("transcript_id", keep="first")


def main():
    parser = argparse.ArgumentParser(description="Attach GSE159044 atlas TSS to CDS annotation by transcript_id.")
    parser.add_argument("--config", default="config/config.yaml")
    parser.add_argument("--cds", default="data/processed/cds_annotation.csv")
    parser.add_argument("--bed", default="data/raw/tss_atlas/GSE159044_eTSS_NCBI_picr.bed.gz")
    parser.add_argument("--meta", default="data/raw/tss_atlas/GSE159044_eTSS_NCBI_picr.meta.tsv.gz")
    parser.add_argument("--output", default="data/processed/tss_annotation.csv")
    parser.add_argument("--all-output", default="data/processed/tss_annotation_all_promoters.csv")
    args = parser.parse_args()

    for p in [args.cds, args.bed, args.meta]:
        if not Path(p).exists():
            raise FileNotFoundError(f"Missing input file: {p}")

    cds = pd.read_csv(args.cds)
    bed = read_tss_bed(Path(args.bed))
    meta = read_tss_meta(Path(args.meta))

    tss = bed.merge(meta, on="tss_id", how="left", suffixes=("", "_meta"))
    merged_all = cds.merge(tss, on="transcript_id", how="inner", suffixes=("_cds", "_tss"))

    if merged_all.empty:
        raise ValueError("No CDS records matched atlas TSS by transcript_id. Check transcript ID parsing.")

    # Use CDS chromosome/strand as authoritative and flag discordant atlas rows.
    merged_all["chromosome_match"] = merged_all["chromosome_cds"] == merged_all["chromosome_tss"]
    merged_all["strand_match"] = merged_all["strand_cds"] == merged_all["strand_tss"]
    merged_all["tss_to_start_codon_distance"] = (merged_all["start_codon_pos"] - merged_all["tss_pos_1based"]).abs()

    best = pick_best_tss(merged_all)

    out_all = Path(args.all_output)
    out = Path(args.output)
    out_all.parent.mkdir(parents=True, exist_ok=True)
    out.parent.mkdir(parents=True, exist_ok=True)
    merged_all.to_csv(out_all, index=False)
    best.to_csv(out, index=False)

    print(f"[OK] All transcript-promoter TSS annotations saved to: {out_all}")
    print(f"[OK] Best TSS annotation saved to: {out}")
    print(f"[INFO] CDS transcripts: {len(cds)}")
    print(f"[INFO] Atlas TSS rows: {len(bed)}")
    print(f"[INFO] Matched transcript-promoter rows: {len(merged_all)}")
    print(f"[INFO] Matched transcripts after best-TSS selection: {len(best)}")
    print(f"[INFO] Chromosome match rate: {merged_all['chromosome_match'].mean():.3f}")
    print(f"[INFO] Strand match rate: {merged_all['strand_match'].mean():.3f}")


if __name__ == "__main__":
    main()
