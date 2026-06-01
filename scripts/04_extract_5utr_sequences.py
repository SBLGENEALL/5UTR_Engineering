import argparse
import gzip
from pathlib import Path

import pandas as pd


COMPLEMENT = str.maketrans("ACGTNacgtn", "TGCANtgcan")


def open_text(path: Path):
    return gzip.open(path, "rt", encoding="utf-8", errors="ignore") if path.suffix == ".gz" else path.open("r", encoding="utf-8", errors="ignore")


def reverse_complement(seq: str) -> str:
    return seq.translate(COMPLEMENT)[::-1].upper()


def read_fasta(path: Path) -> dict[str, str]:
    records = {}
    current_id = None
    chunks = []
    with open_text(path) as handle:
        for line in handle:
            line = line.strip()
            if not line:
                continue
            if line.startswith(">"):
                if current_id is not None:
                    records[current_id] = "".join(chunks).upper()
                current_id = line[1:].split()[0]
                chunks = []
            else:
                chunks.append(line)
    if current_id is not None:
        records[current_id] = "".join(chunks).upper()
    return records


def extract_utr(row, genome: dict[str, str]):
    chrom = row["chromosome_cds"]
    strand = row["strand_cds"]
    tss = int(row["tss_pos_1based"])
    start_codon = int(row["start_codon_pos"])

    if chrom not in genome:
        return None, None, None, "missing_chromosome"

    seq = genome[chrom]

    if strand == "+":
        utr_start = tss
        utr_end = start_codon - 1
        if utr_end < utr_start:
            return utr_start, utr_end, "", "negative_or_zero_length"
        utr_seq = seq[utr_start - 1:utr_end]
    elif strand == "-":
        # On the minus strand, TSS is downstream of the start codon in genomic coordinates.
        utr_start = start_codon + 1
        utr_end = tss
        if utr_end < utr_start:
            return utr_start, utr_end, "", "negative_or_zero_length"
        genomic_seq = seq[utr_start - 1:utr_end]
        utr_seq = reverse_complement(genomic_seq)
    else:
        return None, None, None, "unknown_strand"

    return utr_start, utr_end, utr_seq.upper(), "ok"


def main():
    parser = argparse.ArgumentParser(description="Extract TSS-corrected 5UTR sequences from genome FASTA.")
    parser.add_argument("--config", default="config/config.yaml")
    parser.add_argument("--genome", default="data/reference/GCF_003668045.1_CriGri-PICR_genomic.fna.gz")
    parser.add_argument("--tss", default="data/processed/tss_annotation.csv")
    parser.add_argument("--output", default="data/processed/utr_sequences.csv")
    parser.add_argument("--fasta-output", default="data/processed/utr_sequences.fasta")
    parser.add_argument("--min-length", type=int, default=1)
    parser.add_argument("--max-length", type=int, default=10000)
    args = parser.parse_args()

    genome_path = Path(args.genome)
    tss_path = Path(args.tss)
    if not genome_path.exists():
        raise FileNotFoundError(f"Missing genome FASTA: {genome_path}")
    if not tss_path.exists():
        raise FileNotFoundError(f"Missing TSS annotation: {tss_path}")

    print(f"[INFO] Loading genome FASTA: {genome_path}")
    genome = read_fasta(genome_path)
    print(f"[INFO] Loaded contigs: {len(genome)}")

    tss = pd.read_csv(tss_path)
    rows = []
    for _, row in tss.iterrows():
        utr_start, utr_end, utr_seq, status = extract_utr(row, genome)
        utr_len = len(utr_seq) if isinstance(utr_seq, str) else 0
        out = row.to_dict()
        out.update({
            "utr_start": utr_start,
            "utr_end": utr_end,
            "utr_sequence": utr_seq,
            "utr_length": utr_len,
            "utr_status": status,
        })
        rows.append(out)

    df = pd.DataFrame(rows)
    df["passes_length_filter"] = (
        (df["utr_status"] == "ok")
        & (df["utr_length"] >= args.min_length)
        & (df["utr_length"] <= args.max_length)
    )

    out = Path(args.output)
    fasta_out = Path(args.fasta_output)
    out.parent.mkdir(parents=True, exist_ok=True)
    fasta_out.parent.mkdir(parents=True, exist_ok=True)

    df.to_csv(out, index=False)

    with fasta_out.open("w", encoding="utf-8") as handle:
        for _, r in df[df["passes_length_filter"]].iterrows():
            header = f">{r['transcript_id']}|gene_id={r['gene_id']}|gene={r.get('gene_symbol_cds', '')}|len={r['utr_length']}|strand={r['strand_cds']}"
            handle.write(header + "\n")
            seq = r["utr_sequence"]
            for i in range(0, len(seq), 80):
                handle.write(seq[i:i+80] + "\n")

    print(f"[OK] UTR table saved to: {out}")
    print(f"[OK] UTR FASTA saved to: {fasta_out}")
    print(f"[INFO] Input TSS transcripts: {len(df)}")
    print(f"[INFO] OK UTRs: {(df['utr_status'] == 'ok').sum()}")
    print(f"[INFO] Passing length filter: {df['passes_length_filter'].sum()}")
    print("[INFO] Status counts:")
    print(df["utr_status"].value_counts(dropna=False).to_string())


if __name__ == "__main__":
    main()
