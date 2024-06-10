import pandas as pd
import time

start_time = time.time()

te_list = pd.read_csv('../files/P1_MITE_list.csv', sep='\t').rename(columns={'chr': 'chr_te'})
# te_list = pd.read_csv('../files/P4_MITE_list.csv', sep='\t').rename(columns={'chr': 'chr_te'})

de_and_bins = pd.read_csv('../files/P1_DE_results_comparing_with_genes_in_bins.csv', sep='\t')
# de_and_bins = pd.read_csv('../files/P1_DE_results_comparing_with_genes_in_rest_bins.csv', sep='\t')

# de_and_bins = (pd.read_csv('../files/P1_XLOC_comparing_with_bins_with_genes.csv', sep='\t').
#                rename(columns={'Start': 'start_gene', 'End': 'end_gene'}))

# de_and_bins = (pd.read_csv('../files/P1_XLOC_comparing_with_rest_bins_with_genes.csv', sep='\t').
#                rename(columns={'Start': 'start_gene', 'End': 'end_gene'}))

# de_and_bins = pd.read_csv('../files/P4_DE_results_comparing_with_genes_in_bins.csv', sep='\t')
# de_and_bins = pd.read_csv('../files/P4_DE_results_comparing_with_genes_in_rest_bins.csv', sep='\t')

# de_and_bins = (pd.read_csv('../files/P4_XLOC_comparing_with_bins_with_genes.csv', sep='\t').
#                rename(columns={'Start': 'start_gene', 'End': 'end_gene'}))

# de_and_bins = (pd.read_csv('../files/P4_XLOC_comparing_with_rest_bins_with_genes.csv', sep='\t').
#                rename(columns={'Start': 'start_gene', 'End': 'end_gene'}))

print(te_list.to_string(max_rows=30))
print(de_and_bins.to_string(max_rows=30))

te_in_bins = (de_and_bins.merge(te_list, how='cross').query('(start >= start_bin) & (end <= end_bin) & '
                                                            '(chr == chr_te)').
              drop_duplicates(['start_bin', 'end_bin']).reset_index(drop=True))
print(te_in_bins.to_string(max_rows=30))

te_in_bins_sets = (de_and_bins.merge(te_list, how='cross').query('(start >= start_bin) & (end <= end_bin) & '
                                                                 '(chr == chr_te) & ((te_ins == one_one) & '
                                                                 '(no_te_ins == zero_zero) | (te_ins == zero_zero) & '
                                                                 '(no_te_ins == one_one))').
                   drop_duplicates(['start_bin', 'end_bin', 'chr', 'family', 'start', 'end', 'te_ins', 'no_te_ins']).
                   drop(columns=['start_bin', 'end_bin']).reset_index(drop=True))
print(te_in_bins_sets.to_string())

te_in_genes = de_and_bins.merge(te_list, how='cross').query('(start >= start_gene) & (end <= end_gene) & '
                                                            '(chr == chr_te)').reset_index(drop=True)
print(te_in_genes.to_string(max_rows=30))

te_in_genes_sets = (de_and_bins.merge(te_list, how='cross').
                    query('(start >= start_gene) & (end <= end_gene) & (chr == chr_te) & ((te_ins == one_one) & '
                          '(no_te_ins == zero_zero) | (te_ins == zero_zero) & (no_te_ins == one_one))').
                    reset_index(drop=True))
print(te_in_genes_sets.to_string())

print("--- %s seconds ---" % (time.time() - start_time))
