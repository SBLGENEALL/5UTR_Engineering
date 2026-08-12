# Sort-seq (FACS-seq) 검증 파이프라인

mCherry 리포터(5'UTR variant driving mCherry, GFP-loss로 single-copy 재조합 확인)를
FACS로 MFI 기준 6개 bin + unsorted control로 나눠 gDNA-NGS로 각 variant의 발현을
역추정하고, `07_library_design`에서 예측한 순위(`robust_public_te_rank` 등)와
실측이 맞는지 검증하는 단계입니다. `NGS_LibraryQC`(지난 pDNA amplicon NGS에 쓰인
동일 파이프라인)로 이번에도 raw read를 variant count로 바꾼 뒤,
`5UTR_Engineering`의 새 스크립트(`11_sortseq_bin_expression_score.py`)가 그 결과에
FACS bin 정보(정렬된 세포 수, bin 대표 MFI)를 결합해 variant별 발현 점수를 냅니다.

## 0. 이번 실험이 이전 plasmid NGS run과 다른 점

이전 실행: pDNA(plasmid pool) 1개 샘플만 amplicon-QC → 라이브러리 클로닝 직후
representation(균등하게 잘 들어갔는지)만 확인.

이번 실행: 동일한 reference(같은 `Variant_ID,UTR_sequence` CSV, 같은 primer)로
**8개 샘플**(pDNA 재확인 + unsorted + bin1~bin6)을 한 번에 돌려서, "DNA 단계
representation"과 "세포 내 발현 후 분포"를 같은 좌표계에서 비교합니다.

## 1. Wet-lab에서 반드시 기록해야 하는 것 (스크립트가 요구함)

FACS 소프트웨어(FlowJo/DIVA 등) sort report에서 bin마다:

- `n_cells_sorted`: 실제로 그 bin에 sort된 세포 수 (target 세포 수가 아니라 실측치)
- `bin_mfi_geomean`: 그 gate의 mCherry geometric mean fluorescence (gate % 자체는
  세포수를 모르면 가중치 계산에 못 씀)

이 두 값이 없으면 `11_sortseq_bin_expression_score.py`가 즉시 에러를 내고 멈추도록
만들어 뒀습니다 (조용히 잘못된 값으로 점수를 내는 것을 막기 위함).
`01_pipeline/config/sortseq_bin_design_template.csv`를 복사해서
`01_pipeline/config/sortseq_bin_design.csv`로 채워 넣으세요.

## 2. Step A — amplicon NGS QC (NGS_LibraryQC 저장소, 기존 파이프라인 재사용)

FASTQ 파일명 규칙(Illumina 표준, `<sample>_S#_L00#_R1_001.fastq.gz`)으로
pDNA / unsorted / bin1..bin6 총 8개 샘플을 **같은 input 폴더**에 두고 한 번에 돌립니다.
Reference는 지난 pDNA run과 동일한 `Variant_ID,UTR_sequence` CSV(=
`07_library_design/tables/selected_..._library.csv`에서 뽑은 서열)를 그대로 씁니다.

```bash
cd /data/user/MCET03/03_NGS/01_5UTR_Plasmid  # 또는 이번 run 전용 폴더
python3 amplicon_qc_parallel.py \
  --config config/libraryqc.ini \
  --input-dir raw_data/sortseq_run \
  --outdir results/sortseq_run
```

핵심 산출물: `results/sortseq_run/combined/variant_count_matrix.csv`
(variant × 8 sample raw count matrix). 이걸 그대로 다음 단계 입력으로 씁니다.

## 3. Step B — bin 가중 발현 점수 (5UTR_Engineering 저장소, 신규 스크립트)

```bash
cd /workspace/5utr_engineering  # 프로젝트 루트
cp 01_pipeline/config/sortseq_bin_design_template.csv 01_pipeline/config/sortseq_bin_design.csv
# -> sortseq_bin_design.csv 의 n_cells_sorted / bin_mfi_geomean 를 실측치로 채운 뒤:

python3 01_pipeline/scripts/11_sortseq_bin_expression_score.py \
  --count-matrix /data/.../results/sortseq_run/combined/variant_count_matrix.csv \
  --bin-metadata 01_pipeline/config/sortseq_bin_design.csv \
  --library-design 07_library_design/tables/selected_2000_50_100bp_cluster_diverse_evidence_balanced_library.csv \
  --library-size 2000
```

## 4. 산출물 (`08_sortseq_validation/`)

- `tables/variant_sortseq_expression_score.csv`: variant별 가중 발현 점수
  (`sortseq_log10_mfi_score`), 분포 퍼짐(`sortseq_score_spread`), bin별 추정 세포
  분율, QC flag, (design 테이블을 줬다면) 예측 rank까지 한 줄에.
- `qc/sortseq_qc_summary.txt`: 몇 개가 클로닝 단계에서 빠졌는지 / transfection·selection
  단계에서 빠졌는지 / bin read 부족으로 신뢰 못하는지, pDNA-unsorted 상관, 예측
  rank-실측 점수 상관(Spearman).
- `plots/`: 예측 rank vs 실측 점수 산점도, pDNA vs unsorted 산점도(대표성 유지 확인),
  상/하위 20개 variant의 bin 분포 stacked bar, 전체 점수 분포 히스토그램.

"enrich된 애들이 실제로 발현이 높다"는 주장은 `flag_retained_high_confidence == True`인
variant에서 예측 rank vs `sortseq_log10_mfi_score`의 Spearman r로 정량적으로
뒷받침하십시오. bin 육안 비교만으로 결론 내리지 않는 것이 이 스크립트의 핵심 목적입니다.
