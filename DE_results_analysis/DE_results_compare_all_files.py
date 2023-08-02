import pandas as pd
import ast
import time
from pathlib import Path

start_time = time.time()

genes_in_bins = pd.read_csv('../homozyg_regions/P1_list_genes_in_homozyg_bins.csv', sep='\t')
genes_legend = pd.read_csv('genes.csv', sep=',').rename(columns={'Gene': 'GeneID'})
list_to_de = pd.read_csv('../homozyg_regions/P1_list_to_DE_20more_genes.csv', sep='\t')
list_to_de['df_unique_index'] = list_to_de['df_unique_index'].apply(ast.literal_eval)

final_df = pd.DataFrame()
output_folder = Path("~/Pulpit/burak/DE_results_analysis/P1_DE_results")

for idx in range(0, 2):
    # file_name = f'P1_{idx}_de_result.csv'
    file_name = f'one_one_zero_zero_{idx}_deseq2.xls'
    file_path = output_folder / file_name
    de_results = pd.read_table(file_path)
    # de_results = pd.read_csv(file_path, sep='\t')
    # de_results = pd.read_csv(f'P1_{idx}_de_result.csv', sep='\t')
    print(de_results.to_string(max_rows=30))

    # filter DE results for padj value < 0.05
    de_padj_filter = de_results[de_results['padj'] < 0.05].reset_index(drop=True)
    # print(de_padj_filter.to_string(max_rows=30))

    # comparing set number from file name (<de_results>) with 'P1_list_to_DE_20more_genes.csv' to get <unique_no>
    unique_no = list_to_de.loc[idx, 'df_unique_index']
    print(unique_no, type(unique_no))

    # filter list of genes in bins for appropriate sample set
    genes_in_bins_for_set = genes_in_bins[genes_in_bins['unique_no'].isin(unique_no)]
    print(genes_in_bins_for_set.to_string(max_rows=30))

    # merge filtered DE results with genes legend to identify genes names according to EL10 annotation
    de_padj_filter_merged = de_padj_filter.merge(genes_legend, on=['GeneID']).sort_values(by=['Chromosome', 'Start']). \
        reset_index(drop=True).rename(columns={'Symbol': 'gene_name'})
    # print(de_padj_filter_merged.to_string(max_rows=30))

    # check if there are the same genes in bins for appropriate sample set as identified in DE analysis with padj < 0.05
    de_genes_in_bins = genes_in_bins_for_set.merge(de_padj_filter_merged, on=['gene_name']).drop(
        columns=['genes_count'])
    print(de_genes_in_bins.to_string())

    final_df = final_df.append(de_genes_in_bins)

print(final_df.to_string())

print("--- %s seconds ---" % (time.time() - start_time))





