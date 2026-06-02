# RESULT INTERPRETATION

## 왜 regression보다 classification인가?
현재 label은 대부분 gene-level TE/proteomics 값이고, input은 UTR-level sequence다. 같은 gene의 여러 UTR 후보가 동일 label을 공유할 수 있어 regression scatter에서 같은 true value에 prediction이 세로로 퍼질 수 있다. 따라서 정밀한 순위 예측보다는 high/low 후보군 enrichment로 해석한다.

## protein_residual_rank
`protein_abundance_rank`는 protein이 많이 검출되는 gene이다. 하지만 RNA가 많아서 protein이 많은 gene도 포함될 수 있다. `protein_residual_rank`는 log(protein abundance proxy)를 log(RNA abundance)로 설명한 뒤 residual을 rank화한 값이다. 따라서 RNA 양으로 normalize해도 protein output이 높은 gene을 선별하는 보조 evidence다.

## ROC-AUC
ROC-AUC는 high 후보와 low 후보를 얼마나 잘 구분하는지 보여준다. 무작위 high/low 한 쌍에서 high에 더 높은 score를 줄 확률로 이해하면 된다. 0.5는 랜덤, 0.6 이상은 약한 signal, 0.65–0.70은 proxy label 기준으로 후보 enrichment에 사용할 수 있는 수준이다.

## Jaccard cluster
목적은 동일서열 개수 보고가 아니라, near-duplicate sequence family가 최종 library를 지배하지 않게 하는 것이다. 최종 2000개에서 exact duplicate는 0이어야 하고, 특정 Jaccard cluster에 후보가 과도하게 몰리면 안 된다.

## 최종 2000개 library 해석
`selected_2000_...csv`는 모델 top 2000이 아니다. public TE high, TE classifier high, protein abundance-supported, protein residual-supported, multiomics consensus, exploratory diversity, negative control을 quota로 섞은 실험 검증용 library다.
