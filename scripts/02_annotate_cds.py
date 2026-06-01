import argparse
import gzip
import re
from pathlib import Path
from urllib.parse import unquote

import pandas as pd


PROTEIN_CODING_RNA_TYPES = {"mRNA", "transcript"}


def open_text(path: Path):
    return gzip.open(path, "rt", encoding="utf-8", errors="ignore") if path.suffix == ".gz" else path.open("r", encoding="utf-8", errors="ignore")


def parse_attributes(attr_text: str) -> dict:
    attrs = {}
    for item in attr_text.strip().split(";"):
        if not item:
            continue
        if "=" in item:
            key, value = item.split("=", 1)
            attrs[key] = unquote(value)
    return attrs


def extract_gene_id(attrs: dict) -> str | None:
    dbxref = attrs.get("Dbxref", "")
    match = re.search(r"GeneID:([^,;]+)", dbxref)
    if match:
        return match.group(1)
    return attrs.get("gene") or attrs.get("Name") or attrs.get("ID")


def parse_gff(annotation_path: Path):
    genes = {}
    transcripts = {}
    cds_by_transcript = {}

    with open_text(annotation_path) as handle:
        for line in handle:
            if not line.strip() or line.startswith("#"):
                continue
            parts = line.rstrip("\n").split("\t")
            if len(parts) != 9:
                continue

            chrom, source, feature, start, end, score, strand, phase, attr_text = parts
            start, end = int(start), int(end)
            attrs = parse_attributes(attr_text)

            if feature == "gene":
                gene_record_id = attrs.get("ID")
                gene_id = extract_gene_id(attrs)
                genes[gene_record_id] = {
                    "gene_record_id": gene_record_id,
                    "gene_id": gene_id,
                    "gene_symbol": attrs.get("gene") or attrs.get("Name"),
                    "gene_biotype": attrs.get("gene_biotype"),
                    "chromosome": chrom,
                    "gene_start": start,
                    "gene_end": end,
                    "strand": strand,
                }

            elif feature in PROTEIN_CODING_RNA_TYPES:
                transcript_record_id = attrs.get("ID")
                parent_gene_record_id = attrs.get("Parent")
                gene_info = genes.get(parent_gene_record_id, {})
                transcript_id = attrs.get("transcript_id") or attrs.get("Name") or transcript_record_id
                transcripts[transcript_record_id] = {
                    "transcript_record_id": transcript_record_id,
                    "transcript_id": transcript_id,
                    "parent_gene_record_id": parent_gene_record_id,
                    "gene_id": extract_gene_id(attrs) or gene_info.get("gene_id"),
                    "gene_symbol": attrs.get("gene") or gene_info.get("gene_symbol"),
                    "chromosome": chrom,
                    "transcript_start": start,
                    "transcript_end": end,
                    "strand": strand,
                    "product": attrs.get("product"),
                }

            elif feature == "CDS":
                parent = attrs.get("Parent")
                if not parent:
                    continue
                for parent_id in parent.split(","):
                    cds_by_transcript.setdefault(parent_id, []).append((chrom, start, end, strand, phase, attrs))

    return genes, transcripts, cds_by_transcript


def build_cds_table(transcripts: dict, cds_by_transcript: dict) -> pd.DataFrame:
    rows = []
    for transcript_record_id, cds_list in cds_by_transcript.items():
        tx = transcripts.get(transcript_record_id)
        if tx is None:
            continue
        starts = [x[1] for x in cds_list]
        ends = [x[2] for x in cds_list]
        strand = tx["strand"]
        cds_start = min(starts)
        cds_end = max(ends)
        start_codon_pos = cds_start if strand == "+" else cds_end
        rows.append({
            **tx,
            "cds_start": cds_start,
            "cds_end": cds_end,
            "start_codon_pos": start_codon_pos,
            "cds_exon_count": len(cds_list),
        })
    df = pd.DataFrame(rows)
    if df.empty:
        return df
    df = df.sort_values(["chromosome", "gene_id", "transcript_id"]).reset_index(drop=True)
    return df


def main():
    parser = argparse.ArgumentParser(description="Parse NCBI CriGri-PICR GFF3 and build CDS/start-codon annotation table.")
    parser.add_argument("--config", default="config/config.yaml")
    parser.add_argument("--annotation", default="data/reference/GCF_003668045.1_CriGri-PICR_genomic.gff.gz")
    parser.add_argument("--output", default="data/processed/cds_annotation.csv")
    args = parser.parse_args()

    annotation = Path(args.annotation)
    if not annotation.exists():
        raise FileNotFoundError(f"Missing annotation file: {annotation}")

    _, transcripts, cds_by_transcript = parse_gff(annotation)
    df = build_cds_table(transcripts, cds_by_transcript)

    out = Path(args.output)
    out.parent.mkdir(parents=True, exist_ok=True)
    df.to_csv(out, index=False)

    print(f"[OK] CDS annotation saved to: {out}")
    print(f"[INFO] Protein-coding transcripts with CDS: {len(df)}")
    if not df.empty:
        print(f"[INFO] Unique genes: {df['gene_id'].nunique()}")


if __name__ == "__main__":
    main()
