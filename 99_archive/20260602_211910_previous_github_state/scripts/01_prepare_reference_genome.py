import argparse
import gzip
import json
from pathlib import Path


def exists(path: Path) -> bool:
    return path.exists() and path.is_file()


def count_fasta_records(fasta_path: Path, max_records: int | None = None) -> int:
    opener = gzip.open if fasta_path.suffix == ".gz" else open
    count = 0
    with opener(fasta_path, "rt", encoding="utf-8", errors="ignore") as handle:
        for line in handle:
            if line.startswith(">"):
                count += 1
                if max_records is not None and count >= max_records:
                    break
    return count


def count_annotation_records(annotation_path: Path, max_records: int | None = None) -> int:
    opener = gzip.open if annotation_path.suffix == ".gz" else open
    count = 0
    with opener(annotation_path, "rt", encoding="utf-8", errors="ignore") as handle:
        for line in handle:
            if line.startswith("#") or not line.strip():
                continue
            count += 1
            if max_records is not None and count >= max_records:
                break
    return count


def main():
    parser = argparse.ArgumentParser(description="Validate/register CHO CriGri-PICR reference genome inputs.")
    parser.add_argument("--config", default="config/config.yaml")
    parser.add_argument("--genome", default="data/raw/01_ncbi_genome_annotation/GCF_003668045.1_CriGri-PICR_genomic.fna.gz")
    parser.add_argument("--annotation", default="data/raw/01_ncbi_genome_annotation/GCF_003668045.1_CriGri-PICR_genomic.gff.gz")
    parser.add_argument("--output", default="data/processed/reference_manifest.json")
    args = parser.parse_args()

    genome = Path(args.genome)
    annotation = Path(args.annotation)
    output = Path(args.output)
    output.parent.mkdir(parents=True, exist_ok=True)

    missing = []
    if not exists(genome):
        missing.append(str(genome))
    if not exists(annotation):
        missing.append(str(annotation))
    if missing:
        raise FileNotFoundError("Missing reference input(s):\n" + "\n".join(missing))

    manifest = {
        "genome_fasta": str(genome),
        "annotation_gff_or_gtf": str(annotation),
        "genome_size_bytes": genome.stat().st_size,
        "annotation_size_bytes": annotation.stat().st_size,
        "fasta_record_count_checked": count_fasta_records(genome),
        "annotation_record_count_checked": count_annotation_records(annotation),
        "note": "Reference files validated. Use downstream steps to parse CDS and extract 5UTRs.",
    }

    with output.open("w", encoding="utf-8") as f:
        json.dump(manifest, f, indent=2)

    print(f"[OK] Reference manifest saved to: {output}")
    print(f"[INFO] Genome: {genome}")
    print(f"[INFO] Annotation: {annotation}")


if __name__ == "__main__":
    main()
