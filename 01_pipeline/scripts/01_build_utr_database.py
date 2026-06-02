from pathlib import Path
import json
from collections import defaultdict
import pandas as pd
from common import open_text, parse_attributes, load_fasta_selected, extract_regions_sequence, region_string, write_fasta, ensure_dir

BASE = Path.cwd()
cfg = json.loads((BASE/"01_pipeline/config/project_config.json").read_text(encoding="utf-8"))
GENOME = BASE/cfg["paths"]["genome_fasta"]
ANNOT = BASE/cfg["paths"]["annotation_gff"]
OUT_CSV = BASE/"02_utr_database/tables/ncbi_annotation_5utr_candidates.csv"
OUT_FASTA = BASE/"02_utr_database/fasta/ncbi_annotation_5utr_candidates.fasta"
ensure_dir(OUT_CSV.parent); ensure_dir(OUT_FASTA.parent)

transcripts, exons, cds = {}, defaultdict(list), defaultdict(list)
print("[1] Parse annotation", ANNOT)
with open_text(ANNOT) as f:
    for line in f:
        if not line.strip() or line.startswith("#"): continue
        parts = line.rstrip("\n").split("\t")
        if len(parts) < 9: continue
        seqname, source, ftype, start, end, score, strand, phase, attr = parts[:9]
        start, end = int(start), int(end)
        a = parse_attributes(attr)
        ft = ftype.lower()
        if ft in {"mrna","transcript"}:
            tid = a.get("ID") or a.get("transcript_id") or a.get("Name")
            if not tid: continue
            gid = a.get("Parent") or a.get("gene_id") or a.get("GeneID") or ""
            gname = a.get("gene") or a.get("gene_name") or a.get("Name") or a.get("gene_symbol") or gid or tid
            transcripts[tid] = dict(transcript_id=tid,gene_id=gid,gene_name=gname,seqname=seqname,strand=strand,transcript_start=start,transcript_end=end)
        elif ft == "exon":
            for tid in (a.get("Parent") or a.get("transcript_id") or "").split(","):
                if tid: exons[tid].append((start,end))
        elif ft == "cds":
            for tid in (a.get("Parent") or a.get("transcript_id") or "").split(","):
                if tid: cds[tid].append((start,end))
print("transcripts", len(transcripts), "exon", len(exons), "cds", len(cds))

rows, need = [], set()
for tid, m in transcripts.items():
    if tid not in exons or tid not in cds: continue
    ex = sorted(exons[tid]); cd = sorted(cds[tid]); strand=m["strand"]
    cds_start, cds_end = min(s for s,e in cd), max(e for s,e in cd)
    regs = []
    if strand == "-":
        for s,e in ex:
            us, ue = max(s, cds_end+1), e
            if us <= ue: regs.append((us,ue))
        tss = max(e for s,e in ex)
    else:
        for s,e in ex:
            us, ue = s, min(e, cds_start-1)
            if us <= ue: regs.append((us,ue))
        tss = min(s for s,e in ex)
    if not regs: continue
    need.add(m["seqname"])
    uid = f"{m['gene_name']}_{tid}".replace(" ","_").replace(";","_").replace("/","_")
    rows.append({**m, "utr_id":uid, "annotation_tss":tss, "cds_start":cds_start, "cds_end":cds_end,
                 "utr5_regions":region_string(regs), "utr5_start":min(s for s,e in regs), "utr5_end":max(e for s,e in regs), "tss_confidence":"annotation_only"})
df = pd.DataFrame(rows)
print("[2] candidates before seq extraction", len(df))
genome = load_fasta_selected(GENOME, need)
seqs = []
for _, r in df.iterrows():
    regs = [(int(a), int(b)) for a,b in [x.split("-") for x in str(r.utr5_regions).split(";") if x]]
    seqs.append(extract_regions_sequence(genome, r.seqname, regs, r.strand))
df["utr5_sequence"] = seqs
df = df[df["utr5_sequence"].str.len()>0].copy()
df["utr5_length"] = df["utr5_sequence"].str.len()
df.to_csv(OUT_CSV, index=False)
write_fasta(df, OUT_FASTA, "utr5_sequence", "utr_id")
print("[SAVED]", OUT_CSV, df.shape)
print("[SAVED]", OUT_FASTA)
