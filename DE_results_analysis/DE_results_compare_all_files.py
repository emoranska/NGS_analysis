import pandas as pd
import ast
import time
from pathlib import Path

start_time = time.time()

# genes_in_bins = pd.read_csv('../homozyg_regions/P1_list_genes_in_homozyg_bins.csv', sep='\t')
genes_in_bins = pd.read_csv('../homozyg_regions/P4_list_genes_in_homozyg_bins.csv', sep='\t')
genes_legend = pd.read_csv('genes.csv', sep=',').rename(columns={'Gene': 'GeneID'})

# list_to_de = pd.read_csv('../homozyg_regions/P1_list_to_DE_20more_genes.csv', sep='\t')
list_to_de = pd.read_csv('../homozyg_regions/P4_list_to_DE_20more_genes.csv', sep='\t')
list_to_de['df_unique_index'] = list_to_de['df_unique_index'].apply(ast.literal_eval)

# input_folder = Path("~/Pulpit/Emi_P1")
input_folder = Path("~/Pulpit/Emi_P4")


def de_results_with_genes_in_bins():
    final_df = pd.DataFrame()
    for idx in range(0, 92):
        # read the files with DE results according to path and current sample set number <idx> in loop
        file_name = f'one_one_zero_zero_{idx}_deseq2.xls'
        file_path = input_folder / file_name
        de_results = pd.read_table(file_path)
        # print(de_results.to_string(max_rows=30))

        # filter DE results for padj value < 0.05
        de_padj_filter = de_results[de_results['padj'] < 0.05].reset_index(drop=True)
        # print(de_padj_filter.to_string(max_rows=30))

        # comparing set number from file name (<de_results>) with 'P1_list_to_DE_20more_genes.csv' to get <unique_no>
        unique_no = list_to_de.loc[idx, 'df_unique_index']
        print(unique_no, type(unique_no))

        # filter list of genes in bins for appropriate sample set
        genes_in_bins_for_set = genes_in_bins[genes_in_bins['unique_no'].isin(unique_no)]
        # print(genes_in_bins_for_set.to_string(max_rows=30))

        # merge filtered DE results with genes legend to identify genes names according to EL10 annotation
        de_padj_filter_merged = de_padj_filter.merge(genes_legend, on=['GeneID']).sort_values(by=['Chromosome', 'Start']). \
            reset_index(drop=True).rename(columns={'Symbol': 'gene_name'})
        # print(de_padj_filter_merged.to_string(max_rows=30))

        # check if there are the same genes in bins for appropriate sample set as identified in DE analysis with padj < 0.05
        de_genes_in_bins = genes_in_bins_for_set.merge(de_padj_filter_merged, on=['gene_name']).drop(
            columns=['genes_count', 'Chromosome', 'Strand'])
        de_genes_in_bins['set_no'] = idx
        de_genes_in_bins = de_genes_in_bins.reindex(columns=['set_no', 'unique_no', 'one_one', 'zero_zero', 'chr',
                                                             'start_bin', 'end_bin', 'start_gene', 'end_gene',
                                                             'gene_name', 'gene_strand', 'GeneID', 'Start', 'End',
                                                             'baseMean', 'log2FoldChange', 'lfcSE', 'stat', 'pvalue',
                                                             'padj', 'Biotype', 'Description'])
        # print(de_genes_in_bins.to_string())

        final_df = final_df.append(de_genes_in_bins).reset_index(drop=True)

    print(final_df.to_string())
    # final_df.to_csv('P1_DE_results_comparing_with_genes_in_bins.csv', sep='\t', index=False)
    # final_df.to_csv('P4_DE_results_comparing_with_genes_in_bins.csv', sep='\t', index=False)


de_results_with_genes_in_bins()
print("--- %s seconds ---" % (time.time() - start_time))
