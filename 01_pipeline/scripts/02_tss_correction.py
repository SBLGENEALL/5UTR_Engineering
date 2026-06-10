from pathlib import Path
import argparse
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
OUT_QC = BASE/'03_tss_correction/qc/tss_correction_summary.txt'
ensure_dir(OUT_CSV.parent); ensure_dir(OUT_FASTA.parent); ensure_dir(OUT_QC.parent)


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


def tss_relation(regs, strand, tss, max_extend):
    tss = int(tss)
    if any(s <= tss <= e for s, e in regs):
        return 'inside'
    if not regs or strand not in ['+', '-']:
        return None
    min_start = min(s for s, _ in regs)
    max_end = max(e for _, e in regs)
    if strand == '+':
        dist = min_start - tss
    else:
        dist = tss - max_end
    if 0 < dist <= max_extend:
        return 'upstream'
    return None


def correct_regions(regs, strand, tss, relation):
    if relation == 'inside':
        return trim(regs, strand, tss)
    if relation != 'upstream' or strand not in ['+', '-']:
        return list(regs)
    tss = int(tss)
    out = sorted(regs, key=lambda x: x[0])
    if not out:
        return []
    if strand == '+':
        s, e = out[0]
        out[0] = (min(tss, s), e)
    else:
        s, e = out[-1]
        out[-1] = (s, max(tss, e))
    return out


def choose_tss(cand, regs, strand, max_extend):
    if cand is None:
        return None, None, False
    candidates = []
    wrong_strand_candidate = False
    for _, t in cand.iterrows():
        rel = tss_relation(regs, strand, int(t.tss_pos), max_extend)
        if rel is None:
            continue
        tss_strand = str(t.get('strand', '.'))
        if tss_strand not in [str(strand), '.']:
            wrong_strand_candidate = True
            continue
        score = pd.to_numeric(pd.Series([t.score]), errors='coerce').iloc[0]
        score = float(score) if pd.notna(score) else 1.0
        rel_priority = 0 if rel == 'inside' else 1
        candidates.append((score, -rel_priority, t, rel))
    if not candidates:
        return None, None, wrong_strand_candidate
    candidates.sort(key=lambda x: (x[0], x[1]), reverse=True)
    return candidates[0][2], candidates[0][3], wrong_strand_candidate


ap = argparse.ArgumentParser(description="Correct annotated 5'UTRs with TSS atlas support")
ap.add_argument("--max-extend", type=int, default=300, help="Maximum upstream extension allowed from annotation UTR to atlas TSS")
ap.add_argument("--long-utr-threshold", type=int, default=500, help="QC threshold for unusually long corrected 5'UTRs")
args = ap.parse_args()

df = pd.read_csv(IN_CSV)
rows=[]

if TSS and Path(TSS).exists() and Path(TSS).stat().st_size > 0:
    tss = load_tss(TSS)
    print('UTR rows', len(df), 'TSS rows', len(tss), 'common seqnames', len(set(df.seqname.astype(str)) & set(tss.seqname.astype(str))))
    tss_by_seq = {str(k):v.copy() for k,v in tss.groupby('seqname')}
    for _,r in df.iterrows():
        regs = parse_region_string(r.utr5_regions)
        cand = tss_by_seq.get(str(r.seqname))
        chosen, relation, wrong_strand_candidate = choose_tss(cand, regs, str(r.strand), args.max_extend)
        row = r.to_dict()
        row['tss_qc_seqname_has_tss_atlas'] = cand is not None
        row['tss_qc_chromosome_mismatch'] = cand is None
        row['tss_qc_strand_mismatch'] = bool(wrong_strand_candidate and chosen is None)
        if chosen is not None:
            newregs = correct_regions(regs, str(r.strand), int(chosen.tss_pos), relation)
            if newregs:
                original_regions = region_string(regs)
                corrected_regions = region_string(newregs)
                chosen_score = pd.to_numeric(pd.Series([chosen.score]), errors='coerce').iloc[0]
                row['corrected_tss'] = int(chosen.tss_pos)
                row['tss_signal_score'] = float(chosen_score) if pd.notna(chosen_score) else 1
                row['tss_confidence'] = 'tss_supported_with_signal'
                if corrected_regions == original_regions:
                    row['tss_correction_mode'] = 'unchanged'
                else:
                    row['tss_correction_mode'] = 'trim_to_tss' if relation == 'inside' else 'extend_to_tss'
                row['utr5_regions_tss_corrected'] = corrected_regions
            else:
                row['corrected_tss'] = r.annotation_tss
                row['tss_signal_score'] = np.nan
                row['tss_confidence'] = 'no_tss_match'
                row['tss_correction_mode'] = 'no_tss_match'
                row['utr5_regions_tss_corrected'] = r.utr5_regions
        else:
            row['corrected_tss'] = r.annotation_tss
            row['tss_signal_score'] = np.nan
            row['tss_confidence'] = 'no_tss_match'
            row['tss_correction_mode'] = 'no_tss_match'
            row['utr5_regions_tss_corrected'] = r.utr5_regions
        rows.append(row)
else:
    print('[WARN] TSS atlas missing. Running annotation-only mode: no TSS correction will be applied.')
    for _, r in df.iterrows():
        row = r.to_dict()
        row['corrected_tss'] = r.annotation_tss
        row['tss_signal_score'] = np.nan
        row['tss_confidence'] = 'annotation_only_no_tss_file'
        row['tss_correction_mode'] = 'no_tss_match'
        row['tss_qc_seqname_has_tss_atlas'] = False
        row['tss_qc_chromosome_mismatch'] = True
        row['tss_qc_strand_mismatch'] = False
        row['utr5_regions_tss_corrected'] = r.utr5_regions
        rows.append(row)

out = pd.DataFrame(rows)
genome = load_fasta_selected(GENOME, set(out.seqname.astype(str)))
out['utr5_sequence_tss_corrected'] = [extract_regions_sequence(genome, r.seqname, parse_region_string(r.utr5_regions_tss_corrected), r.strand) for _,r in out.iterrows()]
out['utr5_length_tss_corrected'] = out.utr5_sequence_tss_corrected.str.len()
out['tss_qc_has_N'] = out.utr5_sequence_tss_corrected.fillna('').astype(str).str.contains('N', regex=False)
out['tss_qc_invalid_strand'] = ~out.strand.astype(str).isin(['+', '-'])
out['tss_qc_long_utr'] = out.utr5_length_tss_corrected > args.long_utr_threshold
zero_length_n = int((out.utr5_length_tss_corrected <= 0).sum())
out = out[out.utr5_length_tss_corrected > 0].copy()
out.to_csv(OUT_CSV, index=False)
write_fasta(out, OUT_FASTA, 'utr5_sequence_tss_corrected', 'utr_id')
summary = [
    'TSS correction summary',
    '='*80,
    f'input_rows: {len(df)}',
    f'output_rows: {len(out)}',
    f'zero_length_dropped: {zero_length_n}',
    f'max_extend: {args.max_extend}',
    f'long_utr_threshold: {args.long_utr_threshold}',
    '',
    '[tss_correction_mode]',
    out.tss_correction_mode.value_counts(dropna=False).to_string(),
    '',
    '[tss_confidence]',
    out.tss_confidence.value_counts(dropna=False).to_string(),
    '',
    f'seqname_without_tss_atlas: {int((~out.tss_qc_seqname_has_tss_atlas.fillna(False).astype(bool)).sum())}',
    f'chromosome_mismatch_rows: {int(out.tss_qc_chromosome_mismatch.fillna(False).astype(bool).sum())}',
    f'strand_mismatch_candidates_without_correction: {int(out.tss_qc_strand_mismatch.fillna(False).astype(bool).sum())}',
    f'invalid_strand_rows: {int(out.tss_qc_invalid_strand.fillna(False).astype(bool).sum())}',
    f'N_containing_rows: {int(out.tss_qc_has_N.fillna(False).astype(bool).sum())}',
    f'long_utr_rows: {int(out.tss_qc_long_utr.fillna(False).astype(bool).sum())}',
]
OUT_QC.write_text('\n'.join(summary), encoding='utf-8')
print('[SAVED]', OUT_CSV, out.shape)
print('[SAVED]', OUT_QC)
print(out.tss_confidence.value_counts().to_string())
