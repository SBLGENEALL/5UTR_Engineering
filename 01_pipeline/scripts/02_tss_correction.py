from pathlib import Path
import json
import pandas as pd
import numpy as np
from common import read_table_flexible, load_fasta_selected, parse_region_string, region_string, extract_regions_sequence, write_fasta, ensure_dir

BASE = Path.cwd()
cfg = json.loads((BASE/'01_pipeline/config/project_config.json').read_text(encoding='utf-8'))
GENOME = BASE/cfg['paths']['genome_fasta']
TSS = BASE/cfg['paths'].get('tss_file', '')
IN_CSV = BASE/'02_utr_database/tables/ncbi_annotation_5utr_candidates.csv'
OUT_CSV = BASE/'03_tss_correction/tables/tss_corrected_5utr_database.csv'
OUT_FASTA = BASE/'03_tss_correction/fasta/tss_corrected_5utr_database.fasta'
ensure_dir(OUT_CSV.parent); ensure_dir(OUT_FASTA.parent)


def load_tss(path):
    t = pd.read_csv(path, sep='\t', header=None, comment='#')
    if t.shape[1] < 3:
        raise RuntimeError('TSS bed should have at least 3 columns')
    t = t.iloc[:, :min(6,t.shape[1])].copy()
    t.columns = ['seqname','start0','end','name','score','strand'][:t.shape[1]]
    t['start0'] = pd.to_numeric(t['start0'], errors='coerce')
    t['end'] = pd.to_numeric(t['end'], errors='coerce')
    t['tss_pos'] = ((t['start0']+1+t['end'])/2).round().astype('Int64')
    if 'score' not in t: t['score'] = 1
    if 'strand' not in t: t['strand'] = '.'
    return t[['seqname','tss_pos','strand','score']].dropna(subset=['seqname','tss_pos'])


def trim(regs, strand, tss):
    tss=int(tss); out=[]
    if strand=='-':
        for s,e in regs:
            ns, ne = s, min(e,tss)
            if ns <= ne: out.append((ns,ne))
    else:
        for s,e in regs:
            ns, ne = max(s,tss), e
            if ns <= ne: out.append((ns,ne))
    return out


df = pd.read_csv(IN_CSV)
rows=[]

if TSS and Path(TSS).exists() and Path(TSS).stat().st_size > 0:
    tss = load_tss(TSS)
    print('UTR rows', len(df), 'TSS rows', len(tss), 'common seqnames', len(set(df.seqname.astype(str)) & set(tss.seqname.astype(str))))
    tss_by_seq = {str(k):v.copy() for k,v in tss.groupby('seqname')}
    for _,r in df.iterrows():
        regs = parse_region_string(r.utr5_regions)
        cand = tss_by_seq.get(str(r.seqname))
        chosen = None
        if cand is not None:
            mask = pd.Series(False, index=cand.index)
            for s,e in regs:
                mask |= (cand.tss_pos.astype(float)>=s) & (cand.tss_pos.astype(float)<=e)
            c = cand[mask].copy()
            same = c[(c.strand.astype(str)==str(r.strand)) | (c.strand.astype(str)=='.')]
            if len(same): c = same
            if len(c): chosen = c.sort_values('score', ascending=False).iloc[0]
        row = r.to_dict()
        if chosen is not None:
            newregs = trim(regs, r.strand, int(chosen.tss_pos))
            if newregs:
                row['corrected_tss'] = int(chosen.tss_pos)
                row['tss_signal_score'] = float(chosen.score) if pd.notna(chosen.score) else 1
                row['tss_confidence'] = 'tss_supported_with_signal'
                row['utr5_regions_tss_corrected'] = region_string(newregs)
            else:
                row['corrected_tss'] = r.annotation_tss
                row['tss_signal_score'] = np.nan
                row['tss_confidence'] = 'no_tss_match'
                row['utr5_regions_tss_corrected'] = r.utr5_regions
        else:
            row['corrected_tss'] = r.annotation_tss
            row['tss_signal_score'] = np.nan
            row['tss_confidence'] = 'no_tss_match'
            row['utr5_regions_tss_corrected'] = r.utr5_regions
        rows.append(row)
else:
    print('[WARN] TSS atlas missing. Running annotation-only mode: no TSS correction will be applied.')
    for _, r in df.iterrows():
        row = r.to_dict()
        row['corrected_tss'] = r.annotation_tss
        row['tss_signal_score'] = np.nan
        row['tss_confidence'] = 'annotation_only_no_tss_file'
        row['utr5_regions_tss_corrected'] = r.utr5_regions
        rows.append(row)

out = pd.DataFrame(rows)
genome = load_fasta_selected(GENOME, set(out.seqname.astype(str)))
out['utr5_sequence_tss_corrected'] = [extract_regions_sequence(genome, r.seqname, parse_region_string(r.utr5_regions_tss_corrected), r.strand) for _,r in out.iterrows()]
out['utr5_length_tss_corrected'] = out.utr5_sequence_tss_corrected.str.len()
out = out[out.utr5_length_tss_corrected > 0].copy()
out.to_csv(OUT_CSV, index=False)
write_fasta(out, OUT_FASTA, 'utr5_sequence_tss_corrected', 'utr_id')
print('[SAVED]', OUT_CSV, out.shape)
print(out.tss_confidence.value_counts().to_string())
