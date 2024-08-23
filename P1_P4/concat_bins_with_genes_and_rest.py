import pandas as pd
import numpy as np
import ast

genes_bins = pd.read_csv('../files/P1_list_genes_in_homozyg_bins.csv', sep='\t')
genes_bins_rest = pd.read_csv('../files/P1_list_genes_in_rest_bins_to_add.csv', sep='\t')

# add column with bin length and change its position in df
genes_bins['bin_length'] = genes_bins['end_bin'] - genes_bins['start_bin']
col = genes_bins.pop('bin_length')
genes_bins.insert(3, col.name, col)
print(genes_bins.to_string(max_rows=30))

genes_bins_rest['bin_length'] = genes_bins_rest['end_bin'] - genes_bins_rest['start_bin']
col = genes_bins_rest.pop('bin_length')
genes_bins_rest.insert(3, col.name, col)
print(genes_bins_rest.to_string(max_rows=30))

# create df for all bins with correct bins localisation (according to rest bins - longer) and remove redundant records
all_genes_bins = ((pd.concat([genes_bins, genes_bins_rest], ignore_index=True).sort_values(by=['chr', 'start_bin']).
                  drop_duplicates(subset=['chr', 'start_gene', 'end_gene', 'gene_name', 'unique_no'], keep='last')).
                  reset_index(drop=True))
print(all_genes_bins.to_string(max_rows=100))

genes_count = all_genes_bins.drop_duplicates(subset=['gene_name'])
print(len(genes_count))

sets_list = pd.read_csv('../files/P1_list_to_DE_20more_genes.csv', sep='\t')
print(sets_list.to_string(max_rows=30))


def find_unique_no(sample_set):
    sample_set['one_one'] = sample_set['one_one'].apply(ast.literal_eval)
    sample_set['zero_zero'] = sample_set['zero_zero'].apply(ast.literal_eval)

    # create list of samples to DE with unique_no in separated rows --> to be able to find set_no for every unique_no
    sample_set['df_unique_index'] = sample_set['df_unique_index'].str.split(',')
    sample_set = (sample_set.explode('df_unique_index').rename(columns={'df_unique_index': 'unique_no'}).
                  sort_values(by='unique_no'))
    sample_set['unique_no'] = sample_set['unique_no'].str.replace('[', '').str.replace(']', '')
    sample_set['unique_no'] = sample_set['unique_no'].astype(int)
    print(sample_set.iloc[1, 3], '\n', type(sample_set.iloc[1, 3]), '\n', sample_set.to_string())
    sample_set.to_csv('../files/P1_sets_test.csv', sep='\t', index=False)
    return sample_set


sets_with_unique_no = find_unique_no(sets_list)
all_genes_bins['unique_no'] = all_genes_bins['unique_no'].astype(int)
print(all_genes_bins.iloc[0, 8], type(all_genes_bins.iloc[0, 8]))

# create df for all genes in bins with unique_no for every bin and set_no, one_one and zero_zero for sets typed to DE
all_genes_bins = (all_genes_bins.merge(sets_with_unique_no, how='outer', on='unique_no').
                  drop(columns=['genes_count_y', 'one_one_x', 'zero_zero_x']).
                  rename(columns={'genes_count_x': 'bin_genes_count', 'one_one_y': 'one_one', 'zero_zero_y': 'zero_zero'}).
                  sort_values(by=['chr', 'start_bin', 'start_gene']).
                  reset_index(drop=True))

set_no = all_genes_bins.pop('set_no')
all_genes_bins.insert(8, set_no.name, set_no)
all_genes_bins['set_no'] = all_genes_bins['set_no']

bin_genes_count = all_genes_bins.pop('bin_genes_count')
all_genes_bins.insert(4, bin_genes_count.name, bin_genes_count)
all_genes_bins['bin_genes_count'] = all_genes_bins['bin_genes_count']

print(all_genes_bins.to_string(max_rows=400))

# all_genes_bins = pd.concat([all_genes_bins, sets_list], sort=False).reset_index(drop=True)

# all_genes_bins.to_csv('../files/P1_all_EL10_genes_in_all_bins.csv', sep='\t', index=False)

