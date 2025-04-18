import pandas as pd
import ast
import time
import os
from pathlib import Path
from tpm_calculation import tpm_calc
from P1_P4.concat_bins_with_genes_and_rest import all_genes_in_all_bins

start_time = time.time()

# genes_in_bins = pd.read_csv('../files/P1_list_genes_in_rest_bins_to_add.csv', sep='\t')

genes_in_bins = pd.read_csv('../files/P1_all_EL10_genes_in_all_bins.csv', sep='\t')
# genes_in_bins = all_genes_in_all_bins()
print('Genes in bins', '\n', genes_in_bins.to_string(max_rows=30))

genes_legend = pd.read_csv('../../../Pobrane/genes.csv', sep=',').rename(columns={'Gene': 'GeneID'})

list_to_de = pd.read_csv('../../../Pobrane/P1_list_to_DE_20more_genes.csv', sep='\t')
list_to_de['df_unique_index'] = list_to_de['df_unique_index'].apply(ast.literal_eval)

input_folder = Path("~/Pobrane/Emi_P1")


def de_results_with_genes_in_bins():
    final_df = pd.DataFrame()
    # print(len(os.listdir(input_folder)))

    # set the range according to numer of output files with DE results
    for idx in range(0, 97):
        # read the files with DE results according to path and current sample set number <idx> in loop
        file_name = f'one_one_zero_zero_{idx}_deseq2.xls'
        file_path = input_folder / file_name
        # de_results = pd.read_table(file_path)
        # print(de_results.to_string(max_rows=30))

        # print list with sum of all genes identified in DE for every sample set
        # print(len(de_results))

        # filter DE results for padj value < 0.05
        # de_padj_filter = de_results[de_results['padj'] < 0.05].reset_index(drop=True)
        # print(de_padj_filter.to_string(max_rows=30))

        # print list with sum of genes with padj < 0.05 identified in DE for every sample set
        # print(len(de_padj_filter))

        de_padj_filter = pd.read_table(file_path)

        # merge filtered DE results with genes legend to identify genes names according to EL10 annotation
        de_padj_filter_merged = (de_padj_filter.merge(genes_legend, on=['GeneID']).sort_values
                                 (by=['Chromosome', 'Start']). reset_index(drop=True).rename
                                 (columns={'Symbol': 'gene_name'}))
        # print(de_padj_filter_merged.to_string(max_rows=30))

        # comparing set number from file name (<de_results>) with 'P1_list_to_DE_20more_genes.csv' to get <unique_no>
        unique_no = list_to_de.loc[idx, 'df_unique_index']
        # print(unique_no, type(unique_no))

        # filter list of genes in bins for appropriate sample set
        genes_in_bins_for_set = genes_in_bins[genes_in_bins['unique_no'].isin(unique_no)]
        # print(genes_in_bins_for_set.to_string(max_rows=30))

        # check if there are the same genes in bins for appropriate sample set as identified in DE analysis with
        # padj < 0.05
        # de_genes_in_bins = genes_in_bins_for_set.merge(de_padj_filter_merged, on=['gene_name']).drop(
        #     columns=['Chromosome', 'Strand'])
        # de_genes_in_bins['set_no'] = idx
        # de_genes_in_bins = de_genes_in_bins.reindex(columns=['set_no', 'unique_no', 'one_one', 'zero_zero', 'chr',
        #                                                      'start_bin', 'end_bin', 'start_gene', 'end_gene',
        #                                                      'gene_name', 'gene_strand', 'GeneID', 'Start', 'End',
        #                                                      'baseMean', 'log2FoldChange', 'lfcSE', 'stat', 'pvalue',
        #                                                      'padj', 'Biotype', 'Description'])

        de_genes_in_bins = genes_in_bins_for_set.merge(de_padj_filter_merged, on=['gene_name']).drop(
            columns=['Chromosome', 'Strand'])
        # print(de_genes_in_bins.to_string())

        # pd.append deprecated in pandas > 2.0
        # final_df = final_df.append(de_genes_in_bins).reset_index(drop=True)
        # final_df = final_df.concat(de_genes_in_bins).reset_index(drop=True)
        final_df = pd.concat([final_df, de_genes_in_bins], ignore_index=True)

    # print(final_df.to_string(max_rows=50))
    # final_df.to_csv('P1_DE_results_comparing_with_genes_in_bins.csv', sep='\t', index=False)
    # final_df.to_csv('P4_DE_results_comparing_with_genes_in_bins.csv', sep='\t', index=False)
    # final_df.to_csv('P1_DE_results_comparing_with_genes_in_rest_bins.csv', sep='\t', index=False)
    return final_df


de_df = genes_in_bins.merge(de_results_with_genes_in_bins(), how='outer',
                            on=['chr', 'start_bin', 'end_bin', 'bin_length', 'bin_genes_count', 'start_gene',
                                'end_gene', 'gene_name', 'gene_strand', 'set_no', 'unique_no', 'one_one', 'zero_zero'])

de_df.loc[de_df['log2FoldChange'] > 0, 'de_results'] = 'up'
de_df.loc[de_df['log2FoldChange'] < 0, 'de_results'] = 'down'
de_updown = de_df.pop('de_results')
de_df.insert(14, de_updown.name, de_updown)
de_df['de_results'] = de_df['de_results']

# print(de_df.to_string(max_rows=200))
# de_df.to_csv('../files/P4_all_EL10_genes_in_all_bins_with_DE_results.csv', sep='\t', index=False)

# kallisto = pd.read_csv('../files/kallisto_P4.csv', sep='\t')
# kallisto = pd.read_csv('../files/P1_kallisto_tpm.csv', sep='\t')
kallisto = tpm_calc()

out = de_df.merge(kallisto, on='GeneID', how='left')
# print(out.to_string(max_rows=50))
# print(out.dtypes)

out['one_one'] = out['one_one'].fillna('None')
out['one_one'] = [ast.literal_eval(samples) for samples in out['one_one']]
print(out['one_one'].to_string(max_rows=50))

out['zero_zero'] = out['zero_zero'].fillna('None')
out['zero_zero'] = [ast.literal_eval(samples) for samples in out['zero_zero']]
print(out['zero_zero'].to_string(max_rows=50))

out_one = out.copy()
o = out_one['one_one'].str.join('|').str.get_dummies().astype(bool)
# print(o.to_string(max_rows=50))
out_one.loc[:, o.columns] = out_one.loc[:, o.columns].where(o)
out_one['mean_one_one'] = out_one.loc[:, o.columns].mean(axis=1)
out_one['std_one_one'] = out_one.loc[:, o.columns].std(axis=1)
print(out_one.to_string(max_rows=50))

out_zero = out.copy()
z = out_zero['zero_zero'].str.join('|').str.get_dummies().astype(bool)
# print(z.to_string(max_rows=50))
out_zero.loc[:, z.columns] = out_zero.loc[:, z.columns].where(z)
out_zero['mean_zero_zero'] = out_zero.loc[:, z.columns].mean(axis=1)
out_zero['std_zero_zero'] = out_zero.loc[:, z.columns].std(axis=1)
print(out_zero.to_string(max_rows=50))

k_one_one = out_one.iloc[:, -14:]
k_zero_zero = out_zero.iloc[:, -14:]
k_zero_last2 = k_zero_zero.loc[:, 'mean_zero_zero':]

k_all = k_one_one.fillna(k_zero_zero)
k_all = pd.concat([k_all, k_zero_last2], axis=1)
print(k_all.to_string(max_rows=30))

out_with_kallisto = pd.concat([de_df, k_all], axis=1)
print(out_with_kallisto.to_string(max_rows=200))
# out_with_kallisto.to_csv('../files/P1_all_EL10_genes_in_all_bins_with_DE_and_kallisto_tpm.csv', sep='\t', index=False)


def xloc_checking():
    final_df = pd.DataFrame()

    # set the range according to numer of output files with DE results
    for idx in range(0, 96):
        # read the files with DE results according to path and current sample set number <idx> in loop
        file_name = f'one_one_zero_zero_{idx}_deseq2.xls'
        file_path = input_folder / file_name
        de_results = pd.read_table(file_path)

        # filter DE results for GeneID = XLOC
        de_padj_filter = de_results[de_results.GeneID.str.startswith('XLOC')].reset_index(drop=True)
        # print(de_padj_filter.to_string(max_rows=30))
        print(len(de_padj_filter))

        # merge filtered DE results with genes legend to add columns with coordinates for XLOC genes
        de_padj_filter_merged = de_padj_filter.merge(genes_legend, on=['GeneID']).sort_values(
            by=['Chromosome', 'Start']). \
            reset_index(drop=True).drop(columns=['Symbol', 'Description']).rename(columns={'Chromosome': 'chr'})
        # print(de_padj_filter_merged.to_string(max_rows=30))

        # comparing set number from file name (<de_results>) with 'P1_list_to_DE_20more_genes.csv' to get <unique_no>
        unique_no = list_to_de.loc[idx, 'df_unique_index']
        # print(unique_no, type(unique_no))

        # filter list of bins with genes from EL10 annotation for appropriate sample set
        genes_in_bins_for_set = genes_in_bins[genes_in_bins['unique_no'].isin(unique_no)]
        # print(genes_in_bins_for_set.to_string(max_rows=30))

        # check if in bins for appropriate sample set there are XLOC genes identified in DE analysis
        de_xloc_genes_in_bins = (((genes_in_bins_for_set.merge(de_padj_filter_merged, on=['chr']).query
                                 ('start_bin <= Start <= end_bin')).query('start_bin <= End <= end_bin').
                                 drop(columns=['genes_count', 'start_gene', 'end_gene', 'gene_name', 'gene_strand'])).
                                 drop_duplicates()).reindex(columns=['chr', 'start_bin', 'end_bin', 'set_no',
                                                                     'unique_no', 'one_one', 'zero_zero', 'GeneID',
                                                                     'Start', 'End', 'Strand', 'baseMean',
                                                                     'log2FoldChange', 'lfcSE', 'stat', 'pvalue',
                                                                     'padj', 'Biotype'])
        de_xloc_genes_in_bins['set_no'] = idx
        # print(de_xloc_genes_in_bins.to_string(max_rows=50))

        # pd.append deprecated in pandas > 2.0
        # final_df = final_df.append(de_xloc_genes_in_bins).reset_index(drop=True)
        final_df = pd.concat([final_df, de_xloc_genes_in_bins], ignore_index=True)

        # final_df.to_csv('P1_XLOC_comparing_with_rest_bins_with_genes.csv', sep='\t', index=False)
        # final_df.to_csv('P4_XLOC_comparing_with_bins_with_genes.csv', sep='\t', index=False)

    print(final_df.to_string())


# de_results_with_genes_in_bins()
# xloc_checking()

print("--- %s seconds ---" % (time.time() - start_time))
