from pathlib import Path
import json, gzip

BASE = Path.cwd()
cfg = json.loads((BASE/'01_pipeline/config/project_config.json').read_text(encoding='utf-8'))

print('=== CHO5UTR FINAL NUMBERED PIPELINE INPUT CHECK ===')
ok = True

def check(label, rel, required=True, note=''):
    global ok
    p = BASE / rel
    status = 'OK' if p.exists() and p.stat().st_size > 0 else ('OPTIONAL_MISSING' if not required else 'MISSING')
    print(f'{label:28s} {status:18s} {p} {note}')
    if required and status == 'MISSING':
        ok = False

# Core publicTE inputs
for k, v in cfg['paths'].items():
    required = True  # Final team-release mode: TSS atlas is required.
    note = ''
    check(k, v, required=required, note=note)

# Proteomics minimal TSV can be supplied directly or via included data_assets gz.
prot_tsv = '00_raw_data/05_cho_proteomics/Heffner_minimal.tsv'
prot_csv = '00_raw_data/05_cho_proteomics/Heffner_minimal.csv'
prot_gz_asset = 'data_assets/Heffner_minimal.tsv.gz'
if (BASE/prot_tsv).exists() or (BASE/prot_csv).exists():
    check('Heffner minimal proteomics', prot_tsv, required=False, note='or CSV exists')
elif (BASE/prot_gz_asset).exists():
    print(f"{'Heffner minimal proteomics':28s} {'OK_AS_GZ_ASSET':18s} {BASE/prot_gz_asset} (RUN script will gunzip into raw_data)")
else:
    check('Heffner minimal proteomics', prot_tsv, required=True, note='or data_assets/Heffner_minimal.tsv.gz')
    ok = False

# NCBI gene mapping. gene2accession is essential. gene_info improves symbols and is strongly recommended.
check('NCBI gene2accession', '00_raw_data/05_cho_proteomics/ncbi_gene_mapping/gene2accession.gz', required=True)
check('NCBI gene_info', '00_raw_data/05_cho_proteomics/ncbi_gene_mapping/gene_info.gz', required=False, note='recommended')

if ok:
    print('\nPASS: all required raw-data inputs are present, including TSS atlas.')
else:
    raise SystemExit('\nERROR: required raw-data inputs are missing.')
