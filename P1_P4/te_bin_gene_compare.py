import pandas as pd
import time

start_time = time.time()

te_list = pd.read_csv('../files/P1_MITE_list.csv', sep='\t').rename(columns={'chr': 'chr_te'})
# te_list = pd.read_csv('../files/P4_MITE_list.csv', sep='\t').rename(columns={'chr': 'chr_te'})

de_bins = pd.read_csv('../files/P1_DE_results_comparing_with_genes_in_bins.csv', sep='\t')
de_bins_rest = pd.read_csv('../files/P1_DE_results_comparing_with_genes_in_rest_bins.csv', sep='\t')

# concat DE results for genes in bins and rest bins together
de_and_bins = pd.concat([de_bins, de_bins_rest]).drop_duplicates(subset=['one_one', 'zero_zero', 'chr', 'start_gene',
                                                                         'end_gene', 'gene_name', 'GeneID', 'Start',
                                                                         'End']).reset_index(drop=True)

# count the number of bins with genes identified in DE with padj > 0,5
bins_number = de_and_bins.groupby(['chr', 'start_bin', 'end_bin']).size()
bins_number_chr = bins_number.groupby(['chr']).size()
bins_number_all = bins_number_chr.agg({'': 'sum'})

# count the sum of genes identified in DE with padj > 0,5
genes_in_bins_number = bins_number.agg({'chr': 'sum'})

print(bins_number.to_string())
print(bins_number_chr.to_string())
print(f'Bins number is: {bins_number_all.to_string()}')
print(f'Genes number is: {genes_in_bins_number.to_string()}')
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
print(de_and_bins.to_string(max_rows=60))

# find MITEs in bins
te_in_bins = (de_and_bins.merge(te_list, how='cross').query('(start >= start_bin) & (end <= end_bin) & '
                                                            '(chr == chr_te)').reset_index(drop=True))
print(te_in_bins.to_string(max_rows=30))

# find MITEs in bins with the same sample sets
te_in_bins_sets = (de_and_bins.merge(te_list, how='cross').query('(start >= start_bin) & (end <= end_bin) & '
                                                                 '(chr == chr_te) & ((te_ins == one_one) & '
                                                                 '(no_te_ins == zero_zero) | (te_ins == zero_zero) & '
                                                                 '(no_te_ins == one_one))').
                   drop_duplicates(['start_bin', 'end_bin', 'chr', 'family', 'start', 'end', 'te_ins', 'no_te_ins']).
                   drop(columns=['start_bin', 'end_bin']).reset_index(drop=True))
print(te_in_bins_sets.to_string())

# find MITEs in genes
te_in_genes = de_and_bins.merge(te_list, how='cross').query('(start >= start_gene) & (end <= end_gene) & '
                                                            '(chr == chr_te)').reset_index(drop=True)
print(te_in_genes.to_string(max_rows=30))

# find MITEs in genes with the same sample sets
te_in_genes_sets = (de_and_bins.merge(te_list, how='cross').
                    query('(start >= start_gene) & (end <= end_gene) & (chr == chr_te) & ((te_ins == one_one) & '
                          '(no_te_ins == zero_zero) | (te_ins == zero_zero) & (no_te_ins == one_one))').
                    reset_index(drop=True))
print(te_in_genes_sets.to_string())

print("--- %s seconds ---" % (time.time() - start_time))
