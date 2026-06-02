from __future__ import annotations
from pathlib import Path
import gzip, re
from typing import Dict, Iterable, List, Tuple, Optional
import pandas as pd
import numpy as np


def open_text(path):
    path = str(path)
    if path.endswith('.gz'):
        return gzip.open(path, 'rt', encoding='utf-8', errors='replace')
    return open(path, 'rt', encoding='utf-8', errors='replace')


def read_table_flexible(path: str | Path) -> pd.DataFrame:
    """Read CSV/TSV/gz tables robustly.

    GEO count files are often TSV even when their extension is .txt.gz.
    Try tab before comma and reject one-column parses whose header still
    contains delimiters.
    """
    path = Path(path)
    if path.suffix.lower() in ['.xlsx', '.xls']:
        return pd.read_excel(path, engine='openpyxl')
    attempts = []
    for sep in ['\t', ',', None]:
        for enc in ['utf-8', 'utf-8-sig', 'cp949', 'euc-kr', 'latin1']:
            try:
                df = pd.read_csv(path, sep=sep, engine='python', encoding=enc)
                if df.shape[1] == 1:
                    header = str(df.columns[0])
                    if '\t' in header or ',' in header:
                        attempts.append((sep, enc, 'single-column delimiter header'))
                        continue
                return df
            except Exception as e:
                attempts.append((sep, enc, str(e)[:80]))
    raise RuntimeError(f'Could not read table: {path}; attempts={attempts[:10]}')


def clean_seq(seq) -> str:
    if pd.isna(seq):
        return ''
    return str(seq).upper().replace('U','T').replace(' ','').replace('\n','')


def revcomp(seq: str) -> str:
    return seq.translate(str.maketrans('ACGTNacgtn','TGCANtgcan'))[::-1].upper()


def parse_attributes(attr: str) -> Dict[str, str]:
    attr = str(attr)
    d = {}
    if '=' in attr:
        for part in attr.strip().strip(';').split(';'):
            if '=' in part:
                k,v = part.split('=',1); d[k.strip()] = v.strip()
    else:
        for m in re.finditer(r'(\S+)\s+"([^"]+)"', attr):
            d[m.group(1)] = m.group(2)
    return d


def load_fasta_selected(path: str | Path, wanted_seqnames: Optional[set] = None) -> Dict[str, str]:
    seqs, name, chunks = {}, None, []
    with open_text(path) as f:
        for line in f:
            line = line.rstrip('\n')
            if not line: continue
            if line.startswith('>'):
                if name is not None and (wanted_seqnames is None or name in wanted_seqnames):
                    seqs[name] = ''.join(chunks).upper()
                name = line[1:].split()[0]
                chunks = []
            else:
                if wanted_seqnames is None or name in wanted_seqnames:
                    chunks.append(line.strip())
        if name is not None and (wanted_seqnames is None or name in wanted_seqnames):
            seqs[name] = ''.join(chunks).upper()
    return seqs


def parse_region_string(s: str) -> List[Tuple[int, int]]:
    if pd.isna(s) or not str(s).strip(): return []
    out = []
    for part in str(s).split(';'):
        if part and '-' in part:
            a,b = part.split('-'); out.append((int(a), int(b)))
    return out


def region_string(regions: Iterable[Tuple[int, int]]) -> str:
    return ';'.join(f'{int(a)}-{int(b)}' for a,b in regions if int(a) <= int(b))


def extract_regions_sequence(genome: Dict[str,str], seqname: str, regions: List[Tuple[int,int]], strand: str) -> str:
    if seqname not in genome: return ''
    chrom = genome[seqname]
    parts = []
    ordered = sorted(regions, key=lambda x: x[0], reverse=(strand=='-'))
    for s,e in ordered:
        seg = chrom[s-1:e]
        parts.append(revcomp(seg) if strand=='-' else seg.upper())
    return ''.join(parts).upper()


def write_fasta(df: pd.DataFrame, path: str | Path, seq_col: str, id_col: str = 'utr_id'):
    path = Path(path); path.parent.mkdir(parents=True, exist_ok=True)
    with open(path, 'w', encoding='utf-8') as out:
        for i,r in df.iterrows():
            seq = clean_seq(r.get(seq_col, ''))
            if not seq: continue
            uid = str(r.get(id_col, f'row_{i}')).replace(' ','_').replace('/','_')
            gene = str(r.get('gene_name', 'NA')).replace(' ','_').replace('/','_')
            out.write(f'>{uid}|{gene}|len={len(seq)}\n')
            for j in range(0, len(seq), 80):
                out.write(seq[j:j+80]+'\n')


def ensure_dir(path: str | Path):
    Path(path).mkdir(parents=True, exist_ok=True)


def rank_pct(s):
    return pd.to_numeric(s, errors='coerce').rank(pct=True)


def gc_content(seq):
    seq = clean_seq(seq)
    return (seq.count('G')+seq.count('C'))/len(seq) if len(seq) else np.nan


def find_gene_symbol_col(df: pd.DataFrame) -> str:
    """Find gene symbol/name column in public count tables."""
    candidates = [
        'genesymbol', 'gene_symbol', 'gene symbol', 'GeneSymbol', 'symbol', 'Symbol',
        'gene_name', 'GeneName', 'gene', 'Gene', 'Name'
    ]
    lower = {str(c).lower().replace(' ', '').replace('_',''): c for c in df.columns}
    for cand in candidates:
        key = cand.lower().replace(' ', '').replace('_','')
        if key in lower:
            return lower[key]
    # Fallback: first non-numeric text-like column with gene/symbol in name
    for c in df.columns:
        lc = str(c).lower()
        if 'gene' in lc or 'symbol' in lc:
            return c
    raise KeyError(f'Could not auto-detect gene symbol column. Available={list(df.columns)}')


def find_numeric_cols(df: pd.DataFrame, exclude: List[str] | None = None) -> List[str]:
    exclude = set(exclude or [])
    out = []
    for c in df.columns:
        if c in exclude: continue
        x = pd.to_numeric(df[c], errors='coerce')
        if x.notna().sum() >= max(3, int(0.1*len(df))):
            out.append(c)
    return out


def split_day_cols(cols: List[str]) -> tuple[List[str], List[str]]:
    """Split GSE79512 sample columns into day3/day6.

    Historical project convention:
      RNA day3  = s01-s03
      Ribo day3 = s04-s06
      RNA day6  = s07-s09
      Ribo day6 = s10-s12

    This function returns first three and last three for a given RNA or Ribo table.
    The caller should pass the numeric columns for that table only.
    """
    cols = list(cols)
    if len(cols) >= 6:
        return cols[:3], cols[-3:]
    mid = max(1, len(cols)//2)
    return cols[:mid], cols[mid:]
