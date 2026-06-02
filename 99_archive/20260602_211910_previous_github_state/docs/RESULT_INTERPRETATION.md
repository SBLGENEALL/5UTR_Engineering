# Result Interpretation

## Multi-omics metrics

protein_abundance
Observed protein abundance from proteomics.

protein_residual
Protein abundance remaining after RNA normalization.
Primary ranking metric used for identifying translationally favorable candidates.

protein_te
Protein abundance relative to RNA abundance.

ribo_te
Ribosome occupancy relative to RNA abundance.

## Cluster-aware benchmark

random
Optimistic estimate.

gene_split
Controls gene leakage.

seq_cluster_split
Controls sequence similarity leakage.

gene_seq_cluster_split
Most conservative evaluation and preferred benchmark.

## Final library

The final selected library balances:

- Public TE evidence
- Proteomics evidence
- Protein residual evidence
- Model-supported evidence
- Sequence diversity
- Experimental controls