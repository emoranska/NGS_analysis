import pandas as pd
import numpy as np
import ast

genes_bins = pd.read_csv('../files/P1_list_genes_in_homozyg_bins.csv', sep='\t')
genes_bins_rest = pd.read_csv('../files/P1_list_genes_in_rest_bins_to_add.csv', sep='\t')
sets_list = pd.read_csv('../files/P4_list_to_DE_20more_genes.csv', sep='\t')


def concat_bins_and_rest():
    # add column with bin length and change its position in df
    genes_bins['bin_length'] = genes_bins['end_bin'] - genes_bins['start_bin']
    col = genes_bins.pop('bin_length')
    genes_bins.insert(3, col.name, col)
    # print(genes_bins.to_string(max_rows=30))

    genes_bins_rest['bin_length'] = genes_bins_rest['end_bin'] - genes_bins_rest['start_bin']
    col = genes_bins_rest.pop('bin_length')
    genes_bins_rest.insert(3, col.name, col)
    # print(genes_bins_rest.to_string(max_rows=30))

    # create df for all bins with correct bins localisation (according to rest bins - longer) and remove redundant
    # records
    all_genes_bins = ((pd.concat([genes_bins, genes_bins_rest], ignore_index=True).sort_values(by=['chr', 'start_bin']).
                      drop_duplicates(subset=['chr', 'start_gene', 'end_gene', 'gene_name', 'unique_no'], keep='last')).
                      reset_index(drop=True))
    # print(all_genes_bins.to_string(max_rows=100))
    return all_genes_bins


genes_count = concat_bins_and_rest().drop_duplicates(subset=['gene_name'])
print(len(genes_count))

# to CHECK!!!

# def alternative_to_literal_eval(str_val):
#     if pd.isna(str_val):  # return NaN if value is NaN
#         return np.nan
#     # Remove the square brackets, split on ',' and strip possible
#     # whitespaces between elements
#     # vals = [v.strip() for v in str_val.strip('[]').split(',')]
#     str_val = str_val.replace('[', '')
#     str_val = str_val.replace(']', '')
#
#     # remove duplicates keeping the original order
#     print(list(dict.fromkeys(str_val)))
#     return list(dict.fromkeys(str_val))
#
#     # print(vals)
#     # return list(vals)
#
#
# line = sets_list.loc[1, 'one_one']
# print(type(line), line)
#
# line = line.replace('[', '')
# line = line.replace(']', '')
# print(line)
# # sets_list['one_one'] = sets_list['one_one'].apply(alternative_to_literal_eval)


def find_unique_no(sample_set):
    print('Before:', '\n', sample_set.to_string())
    sample_set['one_one'] = sample_set['one_one'].apply(ast.literal_eval)
    sample_set['zero_zero'] = sample_set['zero_zero'].apply(ast.literal_eval)

    # sample_set['one_one'] = sample_set['one_one'].apply(alternative_to_literal_eval)
    # sample_set['zero_zero'] = sample_set['zero_zero'].apply(alternative_to_literal_eval)
    print('After:', '\n', sample_set.to_string())

    # create list of samples to DE with unique_no in separated rows --> to be able to find set_no for every unique_no
    sample_set['df_unique_index'] = sample_set['df_unique_index'].str.split(',')
    sample_set = (sample_set.explode('df_unique_index').rename(columns={'df_unique_index': 'unique_no'}).
                  sort_values(by='unique_no'))
    sample_set['unique_no'] = sample_set['unique_no'].str.replace('[', '').str.replace(']', '')
    sample_set['unique_no'] = sample_set['unique_no'].astype(int)
    print(sample_set.iloc[1, 3], '\n', type(sample_set.iloc[1, 3]), '\n', sample_set.to_string())
    # sample_set['one_one'] = [repr(x) for x in sample_set['one_one']]
    # sample_set['zero_zero'] = [repr(x) for x in sample_set['zero_zero']]
    # print('Final:', '\n', sample_set.to_string())
    # sample_set.to_csv('../files/P4_sets_test.csv', sep='\t', index=False)
    return sample_set


def all_genes_in_all_bins():
    sets_with_unique_no = find_unique_no(sets_list)
    concat_bins_and_rest()['unique_no'] = concat_bins_and_rest()['unique_no'].astype(int)
    print(concat_bins_and_rest().iloc[0, 8], type(concat_bins_and_rest().iloc[0, 8]))

    # create df for all genes in bins with unique_no for every bin and set_no, one_one and zero_zero for sets to DE
    all_genes_bins = (concat_bins_and_rest().merge(sets_with_unique_no, how='outer', on='unique_no').
                      drop(columns=['genes_count_y', 'one_one_x', 'zero_zero_x']).
                      rename(columns={'genes_count_x': 'bin_genes_count', 'one_one_y': 'one_one',
                                      'zero_zero_y': 'zero_zero'}).
                      sort_values(by=['chr', 'start_bin', 'start_gene']).
                      reset_index(drop=True))

    set_no = all_genes_bins.pop('set_no')
    all_genes_bins.insert(8, set_no.name, set_no)
    all_genes_bins['set_no'] = all_genes_bins['set_no']

    bin_genes_count = all_genes_bins.pop('bin_genes_count')
    all_genes_bins.insert(4, bin_genes_count.name, bin_genes_count)
    all_genes_bins['bin_genes_count'] = all_genes_bins['bin_genes_count']
    print(all_genes_bins.to_string(max_rows=100))

    # all_genes_bins.to_csv('../files/P4_all_EL10_genes_in_all_bins.csv', sep='\t', index=False)
    return all_genes_bins


all_genes_in_all_bins()
