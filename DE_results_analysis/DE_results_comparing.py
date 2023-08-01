import pandas as pd
import ast
import time

genes_in_bins = pd.read_csv('../homozyg_regions/P1_list_genes_in_homozyg_bins.csv', sep='\t')
# de_results = pd.read_csv('P1_0_de_result.csv', sep='\t')

idx = 0
de_results = pd.read_csv(f'P1_{idx}_de_result.csv', sep='\t')

genes_legend = pd.read_csv('genes.csv', sep=',').rename(columns={'Gene': 'GeneID'})

print(genes_in_bins.to_string(max_rows=50))
print(de_results.to_string(max_rows=50))
print(genes_legend.to_string(max_rows=50))

# comparing set number from file name (<de_results>) with 'P1_list_to_DE_20more_genes.csv' to get <unique_no>
list_to_de = pd.read_csv('../homozyg_regions/P1_list_to_DE_20more_genes.csv', sep='\t')
list_to_de['df_unique_index'] = list_to_de['df_unique_index'].apply(ast.literal_eval)
unique_no = list_to_de.loc[idx, 'df_unique_index']
# unique_no = [0, 1]
print(unique_no, type(unique_no))

# filter list of genes in bins for appropriate sample set
genes_in_bins_for_set = genes_in_bins[genes_in_bins['unique_no'].isin(unique_no)]
print(genes_in_bins_for_set.to_string(max_rows=50))

# filter DE results for padj value < 0.05
de_padj_filter = de_results[de_results['padj'] < 0.05].reset_index(drop=True)
print(de_padj_filter.to_string(max_rows=50))

# merge filtered DE results with genes legend to identify genes names according to EL10 annotation
de_padj_filter_merged = de_padj_filter.merge(genes_legend, on=['GeneID']).sort_values(by=['Chromosome', 'Start']).\
    reset_index(drop=True).rename(columns={'Symbol': 'gene_name'})

print(de_padj_filter_merged.to_string(max_rows=50))

# check if there are the same genes in bins for appropriate sample set as identified in DE analysis with padj < 0.05
de_genes_in_bins = genes_in_bins_for_set.merge(de_padj_filter_merged, on=['gene_name']).drop(columns=['genes_count'])
print(de_genes_in_bins.to_string())

'''
What next:
- improve the output to be more clear
- add lines with automatic comparing set number from file name with 'P1_list_to_DE_20more_genes.csv' to get unique_no
- try to find a solution to analyse DE output files for all sample sets for P1 or P4 in loop
'''
