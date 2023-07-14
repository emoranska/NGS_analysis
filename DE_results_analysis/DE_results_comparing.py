import pandas as pd
import time

genes_in_bins = pd.read_csv('../homozyg_regions/P1_list_genes_in_homozyg_bins.csv', sep='\t')
de_results = pd.read_csv('P1_0_de_result.csv', sep='\t')
genes_legend = pd.read_csv('genes.csv', sep=',').rename(columns={'Gene': 'GeneID'})

print(genes_in_bins.to_string(max_rows=50))
print(de_results.to_string(max_rows=50))
print(genes_legend.to_string(max_rows=50))

unique_no = [0, 1]
genes_in_bins_for_set = genes_in_bins[genes_in_bins['unique_no'].isin(unique_no)]
print(genes_in_bins_for_set.to_string(max_rows=50))

de_padj_filter = de_results[de_results['padj'] < 0.05].reset_index(drop=True)
print(de_padj_filter.to_string(max_rows=50))

de_padj_filter_merged = de_padj_filter.merge(genes_legend, on=['GeneID']).sort_values(by=['Chromosome', 'Start']).\
    reset_index(drop=True).rename(columns={'Symbol': 'gene_name'})

print(de_padj_filter_merged.to_string(max_rows=50))

de_genes_in_bins = genes_in_bins_for_set.merge(de_padj_filter_merged, on=['gene_name'])
print(de_genes_in_bins.to_string())
