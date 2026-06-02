from pathlib import Path
import json
import pandas as pd
import numpy as np
from common import read_table_flexible, rank_pct, gc_content, ensure_dir, find_gene_symbol_col, find_numeric_cols

BASE = Path.cwd()
cfg = json.loads((BASE/'01_pipeline/config/project_config.json').read_text(encoding='utf-8'))
UTR = BASE/'03_tss_correction/tables/tss_corrected_5utr_database.csv'
RNA = BASE/cfg['paths']['rna_count_table']
RIBO = BASE/cfg['paths']['ribo_count_table']
OUT = BASE/'04_te_labeling/tables/tss_corrected_5utr_robust_public_te_labels.csv'
READY = BASE/'04_te_labeling/tables/tss_corrected_5utr_50_100bp_training_ready.csv'
QC = BASE/'04_te_labeling/qc/robust_public_te_mapping_summary.txt'
ensure_dir(OUT.parent); ensure_dir(QC.parent)


def norm_symbol(x): return '' if pd.isna(x) else str(x).strip().upper()

def as_list(x):
    if x is None: return []
    if isinstance(x, list): return x
    if isinstance(x, tuple): return list(x)
    if isinstance(x, str):
        if x.lower() == 'auto': return 'auto'
        if ',' in x: return [v.strip() for v in x.split(',') if v.strip()]
        return [x]
    return [x]


def resolve_day_cols(df, gene_col, table_kind, c_day3, c_day6):
    d3 = as_list(c_day3); d6 = as_list(c_day6)
    if d3 != 'auto' and d6 != 'auto':
        return d3, d6
    nums = find_numeric_cols(df, exclude=[gene_col])
    # Project-specific GSE79512 convention. If the canonical sample names exist, use them.
    if table_kind == 'RNA':
        auto3 = [c for c in ['s01','s02','s03'] if c in df.columns]
        auto6 = [c for c in ['s07','s08','s09'] if c in df.columns]
    else:
        auto3 = [c for c in ['s04','s05','s06'] if c in df.columns]
        auto6 = [c for c in ['s10','s11','s12'] if c in df.columns]
    if len(auto3) == 3 and len(auto6) == 3:
        guess3, guess6 = auto3, auto6
    elif len(nums) >= 6:
        guess3, guess6 = nums[:3], nums[-3:]
    else:
        mid = max(1, len(nums)//2)
        guess3, guess6 = nums[:mid], nums[mid:]
    if d3 != 'auto': guess3 = d3
    if d6 != 'auto': guess6 = d6
    print(f'[AUTO {table_kind}] gene_col={gene_col}; day3={guess3}; day6={guess6}')
    return as_list(guess3), as_list(guess6)


def check_columns(df, cols, label):
    missing = [c for c in cols if c not in df.columns]
    if missing:
        raise KeyError(f'{label} missing columns: {missing}\nAvailable columns: {list(df.columns)}')


def cpm(df, cols):
    out=df.copy()
    for c in cols:
        out[c]=pd.to_numeric(out[c], errors='coerce').fillna(0)
        lib=out[c].sum()
        out[c+'_cpm']=out[c]/lib*1e6 if lib>0 else 0
    return out


def residual_rank(log_rna, log_ribo):
    valid=np.isfinite(log_rna)&np.isfinite(log_ribo)
    res=np.full(len(log_rna), np.nan)
    if valid.sum()>10:
        slope, intercept=np.polyfit(log_rna[valid], log_ribo[valid], 1)
        res[valid]=log_ribo[valid]-(slope*log_rna[valid]+intercept)
    return pd.Series(res).rank(pct=True).values

cols=cfg['columns']; cf=cfg['candidate_filter']; tf=cfg['training_filter']; w=cfg['robust_public_te_score']
utr=pd.read_csv(UTR)
utr['gene_symbol_key']=utr['gene_name'].apply(norm_symbol)
rna=read_table_flexible(RNA); ribo=read_table_flexible(RIBO)
rna_gene = find_gene_symbol_col(rna) if cols.get('rna_gene_symbol') == 'auto' else cols['rna_gene_symbol']
ribo_gene = find_gene_symbol_col(ribo) if cols.get('ribo_gene_symbol') == 'auto' else cols['ribo_gene_symbol']
rna_day3, rna_day6 = resolve_day_cols(rna, rna_gene, 'RNA', cols.get('rna_day3','auto'), cols.get('rna_day6','auto'))
ribo_day3, ribo_day6 = resolve_day_cols(ribo, ribo_gene, 'RIBO', cols.get('ribo_day3','auto'), cols.get('ribo_day6','auto'))
check_columns(rna, [rna_gene]+rna_day3+rna_day6, 'RNA')
check_columns(ribo, [ribo_gene]+ribo_day3+ribo_day6, 'Ribo')

rna['gene_symbol_key']=rna[rna_gene].apply(norm_symbol); ribo['gene_symbol_key']=ribo[ribo_gene].apply(norm_symbol)
rna=cpm(rna, rna_day3+rna_day6); ribo=cpm(ribo, ribo_day3+ribo_day6)
rs=pd.DataFrame({
    'gene_symbol_key':rna.gene_symbol_key,
    'rna_day3_raw_mean':rna[rna_day3].apply(pd.to_numeric, errors='coerce').mean(axis=1),
    'rna_day6_raw_mean':rna[rna_day6].apply(pd.to_numeric, errors='coerce').mean(axis=1),
    'rna_day3_cpm_mean':rna[[c+'_cpm' for c in rna_day3]].mean(axis=1),
    'rna_day6_cpm_mean':rna[[c+'_cpm' for c in rna_day6]].mean(axis=1),
}).groupby('gene_symbol_key', as_index=False).mean(numeric_only=True)
bs=pd.DataFrame({
    'gene_symbol_key':ribo.gene_symbol_key,
    'ribo_day3_raw_mean':ribo[ribo_day3].apply(pd.to_numeric, errors='coerce').mean(axis=1),
    'ribo_day6_raw_mean':ribo[ribo_day6].apply(pd.to_numeric, errors='coerce').mean(axis=1),
    'ribo_day3_cpm_mean':ribo[[c+'_cpm' for c in ribo_day3]].mean(axis=1),
    'ribo_day6_cpm_mean':ribo[[c+'_cpm' for c in ribo_day6]].mean(axis=1),
}).groupby('gene_symbol_key', as_index=False).mean(numeric_only=True)
df=utr.merge(rs,on='gene_symbol_key',how='left').merge(bs,on='gene_symbol_key',how='left')
pc=1e-6
df['te_day3_norm']=(df.ribo_day3_cpm_mean+pc)/(df.rna_day3_cpm_mean+pc)
df['te_day6_norm']=(df.ribo_day6_cpm_mean+pc)/(df.rna_day6_cpm_mean+pc)
df['te_mean_norm']=df[['te_day3_norm','te_day6_norm']].mean(axis=1)
df['log2_te_day3_norm']=np.log2(df.te_day3_norm+pc)
df['log2_te_day6_norm']=np.log2(df.te_day6_norm+pc)
df['log2_te_mean_norm']=np.log2(df.te_mean_norm+pc)
df['mean_TE_rank']=rank_pct(df.log2_te_mean_norm)
df['day3_TE_rank']=rank_pct(df.log2_te_day3_norm)
df['day6_TE_rank']=rank_pct(df.log2_te_day6_norm)
df['day_consensus_TE_rank']=df[['day3_TE_rank','day6_TE_rank']].min(axis=1)
df['ribo_abundance_rank']=rank_pct(np.log2((df.ribo_day3_cpm_mean+df.ribo_day6_cpm_mean)/2+pc))
avg_rna=((df.rna_day3_cpm_mean+df.rna_day6_cpm_mean)/2).fillna(0).values
avg_ribo=((df.ribo_day3_cpm_mean+df.ribo_day6_cpm_mean)/2).fillna(0).values
df['residual_TE_rank']=residual_rank(np.log2(avg_rna+pc), np.log2(avg_ribo+pc))
if 'tss_confidence' in df.columns:
    df['tss_confidence_score']=(df.tss_confidence.astype(str)=='tss_supported_with_signal').astype(float)
else:
    df['tss_confidence_score']=0.0
df['robust_public_te_score']=(w['mean_te_rank_weight']*df.mean_TE_rank.fillna(0)+w['day_consensus_rank_weight']*df.day_consensus_TE_rank.fillna(0)+w['residual_te_rank_weight']*pd.Series(df.residual_TE_rank).fillna(0)+w['ribo_abundance_rank_weight']*df.ribo_abundance_rank.fillna(0)+w['tss_confidence_weight']*df.tss_confidence_score.fillna(0))
df['robust_public_te_rank']=rank_pct(df.robust_public_te_score)
seq=df.utr5_sequence_tss_corrected.fillna('').astype(str)
df['utr5_length_final']=seq.str.len()
df['gc_content']=seq.apply(gc_content)
df['uaug_count']=seq.str.count('ATG')
df['has_rna_label']=df.rna_day3_raw_mean.notna() & df.rna_day6_raw_mean.notna()
df['has_ribo_label']=df.ribo_day3_raw_mean.notna() & df.ribo_day6_raw_mean.notna()
df['primary_length_50_100']=df.utr5_length_final.between(cf['primary_min_length'], cf['primary_max_length'])
df['train_length_40_200']=df.utr5_length_final.between(40,200)
df['gc_pass']=df.gc_content.between(cf['gc_min'], cf['gc_max'])
df['uaug_pass']=df.uaug_count <= int(cf['allow_uaug_max'])
ready=(df.has_rna_label & df.has_ribo_label & df.primary_length_50_100 & df.gc_pass & df.uaug_pass & (df.rna_day3_raw_mean>=tf['rna_day3_raw_mean_min']) & (df.rna_day6_raw_mean>=tf['rna_day6_raw_mean_min']) & (df.ribo_day3_raw_mean>=tf['ribo_day3_raw_mean_min']) & (df.ribo_day6_raw_mean>=tf['ribo_day6_raw_mean_min']) & df.robust_public_te_rank.notna())
df['training_ready_50_100bp']=ready
df.to_csv(OUT,index=False)
df[df.training_ready_50_100bp].to_csv(READY,index=False)
summary=[
    'Robust public TE mapping summary', '='*80,
    f'rna_gene_col: {rna_gene}', f'ribo_gene_col: {ribo_gene}',
    f'rna_day3: {rna_day3}', f'rna_day6: {rna_day6}',
    f'ribo_day3: {ribo_day3}', f'ribo_day6: {ribo_day6}',
    f'UTR_rows: {len(df)}', f'training_ready_50_100bp: {int(ready.sum())}', ''
]
for c in ['has_rna_label','has_ribo_label','primary_length_50_100','train_length_40_200','gc_pass','uaug_pass','training_ready_50_100bp']:
    summary.append('\n'+c+'\n'+df[c].value_counts(dropna=False).to_string())
QC.write_text('\n'.join(summary), encoding='utf-8')
print('[SAVED]', OUT, df.shape)
print('[SAVED]', READY, int(ready.sum()))
print('[SAVED]', QC)
