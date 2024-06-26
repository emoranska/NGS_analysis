import pandas as pd
import time

start_time = time.time()

mite_list_p1 = pd.read_csv('../files/P1_MITE_list.csv', sep='\t').rename(columns={'chr': 'chr_te'})
mite_list_p4 = pd.read_csv('../files/P4_MITE_list.csv', sep='\t').rename(columns={'chr': 'chr_te'})

# print(mite_list_p1.to_string(max_rows=30))
print(mite_list_p4.to_string(max_rows=30))

de_bins_p1 = pd.read_csv('../files/P1_DE_results_comparing_with_genes_in_bins.csv', sep='\t')
de_bins_rest_p1 = pd.read_csv('../files/P1_DE_results_comparing_with_genes_in_rest_bins.csv', sep='\t')

de_bins_p4 = pd.read_csv('../files/P4_DE_results_comparing_with_genes_in_bins.csv', sep='\t')
de_bins_rest_p4 = pd.read_csv('../files/P4_DE_results_comparing_with_genes_in_rest_bins.csv', sep='\t')

xloc_bins_p1 = pd.read_csv('../files/P1_XLOC_comparing_with_bins_with_genes.csv', sep='\t')


xloc_rest_bins_p1 = pd.read_csv('../files/P1_XLOC_comparing_with_rest_bins_with_genes.csv', sep='\t')


xloc_bins_p4 = pd.read_csv('../files/P4_XLOC_comparing_with_bins_with_genes.csv', sep='\t')

xloc_rest_bins_p4 = pd.read_csv('../files/P4_XLOC_comparing_with_rest_bins_with_genes.csv', sep='\t')


def concat_and_count_bins(bins, rest_bins):
    # concat df with bins compared with DE results and rest bins into one
    bins_and_rest = pd.concat([bins, rest_bins]).drop_duplicates(
        subset=['one_one', 'zero_zero', 'chr', 'GeneID', 'Start',
                'End']).reset_index(drop=True)

    # count the number of bins with genes identified in DE with padj > 0,5
    bins_number = bins_and_rest.groupby(['chr', 'start_bin', 'end_bin']).size()
    bins_number_chr = bins_number.groupby(['chr']).size()
    bins_number_all = bins_number_chr.agg({'': 'sum'})

    # count the sum of genes identified in DE with padj > 0,5
    genes_in_bins_number = bins_number.agg({'chr': 'sum'})

    print(bins_number.to_string())
    print(bins_number_chr.to_string())
    print(f'Bins number is: {bins_number_all.iloc[0]}')
    print(f'Genes number is: {genes_in_bins_number.iloc[0]}')
    return bins_and_rest


# de_and_bins_p1 = concat_and_count_bins(de_bins_p1, de_bins_rest_p1)
# xloc_and_bins_p1 = concat_and_count_bins(xloc_bins_p1, xloc_rest_bins_p1)

# de_and_bins_p4 = concat_and_count_bins(de_bins_p4, de_bins_rest_p4)
xloc_and_bins_p4 = concat_and_count_bins(xloc_bins_p4, xloc_rest_bins_p4)


def te_in_bins_and_genes(te_list, genes_and_bins):
    # find MITEs in bins
    te_in_bins = (genes_and_bins.merge(te_list, how='cross').query('(start >= start_bin) & (end <= end_bin) '
                                                                   '& (chr == chr_te)').reset_index(drop=True))
    print(te_in_bins.to_string(max_rows=30))

    # find MITEs in bins with the same sample sets
    te_in_bins_sets = (genes_and_bins.merge(te_list, how='cross').
                       query('(start >= start_bin) & (end <= end_bin) & (chr == chr_te) & ((te_ins == one_one) & '
                             '(no_te_ins == zero_zero) | (te_ins == zero_zero) & (no_te_ins == one_one))').
                       drop_duplicates(['start_bin', 'end_bin', 'chr', 'family', 'start', 'end', 'te_ins', 'no_te_ins'])
                       .drop(columns=['start_bin', 'end_bin']).reset_index(drop=True))
    print(te_in_bins_sets.to_string())

    # find MITEs in genes
    te_in_genes = (genes_and_bins.merge(te_list, how='cross').
                   query('(start >= Start) & (end <= End) & (chr == chr_te)').reset_index(drop=True))
    print(te_in_genes.to_string(max_rows=30))

    # find MITEs in genes with the same sample sets
    te_in_genes_sets = (genes_and_bins.merge(te_list, how='cross').
                        query('(start >= Start) & (end <= End) & (chr == chr_te) & ((te_ins == one_one) & '
                              '(no_te_ins == zero_zero) | (te_ins == zero_zero) & (no_te_ins == one_one))').
                        reset_index(drop=True))

    print(te_in_genes_sets.to_string())
    print(f'Sum of MITEs with 3 and 3 samples in set: {len(te_list)}')
    print(f'Sum of MITEs in bins: {len(te_in_bins)} and with the same sample sets: {len(te_in_bins_sets)}')
    print(f'Sum of MITEs in genes: {len(te_in_genes)} and with the same sample sets: {len(te_in_genes_sets)}')


# te_in_bins_and_genes(mite_list_p1, de_and_bins_p1)
# te_in_bins_and_genes(mite_list_p1, xloc_and_bins_p1)
# te_in_bins_and_genes(mite_list_p4, de_and_bins_p4)
te_in_bins_and_genes(mite_list_p4, xloc_and_bins_p4)

print("--- %s seconds ---" % (time.time() - start_time))
