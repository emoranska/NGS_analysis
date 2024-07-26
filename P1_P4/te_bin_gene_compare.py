import pandas as pd
import time
import ast

start_time = time.time()

mite_list_p1 = pd.read_csv('../files/P1_MITEs_with_bins_sets.csv', sep='\t').rename(columns={'chr': 'chr_te'})
mite_list_p4 = pd.read_csv('../files/P4_MITEs_with_bins_sets_2.csv', sep='\t').rename(columns={'chr': 'chr_te'})

# print(mite_list_p1.to_string(max_rows=30))
# print(mite_list_p4.to_string(max_rows=30))

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
    bins_and_rest = pd.concat([bins, rest_bins]).sort_values(by=['chr', 'start_bin', 'end_bin']).drop_duplicates(
        subset=['one_one', 'zero_zero', 'chr', 'start_bin', 'end_bin', 'GeneID', 'Start',
                'End']).reset_index(drop=True)

    bins_and_rest['bin_length'] = bins_and_rest['end_bin'] - bins_and_rest['start_bin']
    col = bins_and_rest.pop('bin_length')
    bins_and_rest.insert(7, col.name, col)
    print(bins_and_rest.to_string(max_rows=150))

    bins_and_rest = (bins_and_rest.sort_values(by=['chr', 'start_bin', 'end_bin', 'bin_length']).
                     drop_duplicates(subset=['one_one', 'zero_zero', 'chr', 'Start', 'End', 'GeneID'],
                                     keep='last').reset_index(drop=True))
    print(bins_and_rest.to_string(max_rows=150))

    # count the number of bins with genes identified in DE with padj < 0.05
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


def clean_sample_list(col):
    col = col.apply(ast.literal_eval)
    col = [', '.join(map(str, l)) for l in col]
    return col


def te_in_bins_and_genes(te_list, genes_and_bins):

    # # find MITEs in bins with the same sample sets for annotated genes
    # te_in_bins_sets = ((genes_and_bins.merge(te_list, how='cross').
    #                    query('(start >= start_bin) & (end <= end_bin) & (chr == chr_te) & ((te_ins == one_one) & '
    #                          '(no_te_ins == zero_zero) | (te_ins == zero_zero) & (no_te_ins == one_one))').
    #                    drop_duplicates(['start_bin', 'end_bin', 'chr', 'family', 'start', 'end', 'te_ins', 'no_te_ins'])
    #                    .drop(columns=['start_gene', 'end_gene', 'GeneID', 'baseMean', 'stat', 'pvalue', 'padj', 'Unnamed: 0', 'chr_te']).
    #                    rename(columns={'chr': 'chromosome', 'Start': 'start_gene', 'End': 'end_gene',
    #                                    'gene_name': 'gene_symbol', 'Strand': 'gene_strand', 'Biotype': 'biotype',
    #                                    'Description': 'gene_function', 'start': 'start_te', 'end': 'end_te'})).
    #                    sort_values(by=['chromosome', 'start_te']).reset_index(drop=True))
    #
    # te_in_bins_sets = te_in_bins_sets[['chromosome', 'family', 'start_te', 'end_te', 'te_ins', 'no_te_ins', 'set_no',
    #                                    'unique_no', 'start_bin', 'end_bin',
    #                                    'gene_symbol', 'start_gene', 'end_gene', 'gene_strand', 'one_one', 'zero_zero',
    #                                    'log2FoldChange', 'lfcSE', 'biotype', 'gene_function'
    #                                    ]]

    # find MITEs in bins with the same sample sets for XLOC
    te_in_bins_sets = ((genes_and_bins.merge(te_list, how='cross').
                       query('(start >= start_bin) & (end <= end_bin) & (chr == chr_te) & ((te_ins == one_one) & '
                             '(no_te_ins == zero_zero) | (te_ins == zero_zero) & (no_te_ins == one_one))').
                       drop_duplicates(['start_bin', 'end_bin', 'chr', 'family', 'start', 'end', 'te_ins', 'no_te_ins'])
                       .drop(columns=['baseMean', 'stat', 'pvalue', 'padj', 'Unnamed: 0', 'chr_te']).
                       rename(columns={'chr': 'chromosome', 'Start': 'start_gene', 'End': 'end_gene',
                                       'GeneID': 'gene_symbol', 'Strand': 'gene_strand', 'Biotype': 'biotype',
                                       'start': 'start_te', 'end': 'end_te'})).
                       sort_values(by=['chromosome', 'start_te']).reset_index(drop=True))

    te_in_bins_sets = te_in_bins_sets[['chromosome', 'family', 'start_te', 'end_te', 'te_ins', 'no_te_ins', 'set_no',
                                       'unique_no', 'start_bin', 'end_bin',
                                       'gene_symbol', 'start_gene', 'end_gene', 'gene_strand', 'one_one', 'zero_zero',
                                       'log2FoldChange', 'lfcSE', 'biotype'
                                       ]]

    te_in_bins_sets['te_ins'] = clean_sample_list(te_in_bins_sets['te_ins'])
    te_in_bins_sets['no_te_ins'] = clean_sample_list(te_in_bins_sets['no_te_ins'])
    te_in_bins_sets['one_one'] = clean_sample_list(te_in_bins_sets['one_one'])
    te_in_bins_sets['zero_zero'] = clean_sample_list(te_in_bins_sets['zero_zero'])

    print(te_in_bins_sets.to_string(max_rows=50))
    unique_bins_count = te_in_bins_sets.drop_duplicates(subset=['start_bin', 'end_bin'])
    print(f'Number of bins with MITEs: {len(unique_bins_count)}')

    # te_in_bins_sets.to_csv('../files/P4_MITEs_in_bins_1.csv', sep='\t')

    # find MITEs in genes with the same sample sets
    te_in_genes_sets = (genes_and_bins.merge(te_list, how='cross').
                        query('(start >= Start) & (end <= End) & (chr == chr_te) & ((te_ins == one_one) & '
                              '(no_te_ins == zero_zero) | (te_ins == zero_zero) & (no_te_ins == one_one))').
                        reset_index(drop=True))

    # # find MITEs upstream/downstream 2000 bp from genes with the same sample sets
    # te_in_genes_sets = (genes_and_bins.merge(te_list, how='cross').
    #                     query('(((start >= Start - 2000) & (end < Start) & (chr == chr_te)) | '
    #                           '((end <= End + 2000) & (start > End) & (chr == chr_te))) & ((te_ins == one_one) & '
    #                           '(no_te_ins == zero_zero) | (te_ins == zero_zero) & (no_te_ins == one_one))').
    #                     reset_index(drop=True))

    print(f'Sum of MITEs with the same 3 and 3 samples in set as sets to DE (from bins): {len(te_list)}')
    print(f'Sum of MITEs in bins with the same sample sets: {len(te_in_bins_sets)}')

    unique_bins_count_m = te_in_genes_sets.drop_duplicates(subset=['start_bin', 'end_bin'])
    print(f'Number of bins with MITEs in DE genes: {len(unique_bins_count_m)}')

    # # cleaning the output for annotated genes
    # te_in_genes_sets = ((te_in_genes_sets.drop(
    #     columns=['start_gene', 'end_gene', 'GeneID', 'baseMean', 'stat', 'pvalue', 'padj', 'chr_te']).
    #                     rename(columns={'chr': 'chromosome', 'Start': 'start_gene', 'End': 'end_gene',
    #                                     'gene_name': 'gene_symbol', 'Strand': 'gene_strand', 'Biotype': 'biotype',
    #                                     'Description': 'gene_function', 'start': 'start_te',
    #                                     'end': 'end_te'})).sort_values(by=['chromosome', 'start_te']).
    #                     reset_index(drop=True))
    #
    # te_in_genes_sets = te_in_genes_sets[['chromosome', 'family', 'start_te', 'end_te', 'te_ins', 'no_te_ins',
    #                                      'gene_symbol', 'start_gene', 'end_gene', 'gene_strand', 'one_one', 'zero_zero',
    #                                      'log2FoldChange', 'lfcSE', 'biotype', 'gene_function', 'start_bin', 'end_bin',
    #                                      'set_no', 'unique_no']]

    # cleaning the output for XLOC genes
    te_in_genes_sets = ((te_in_genes_sets.drop(
        columns=['baseMean', 'stat', 'pvalue', 'padj', 'chr_te']).
                        rename(columns={'chr': 'chromosome', 'GeneID': 'gene_symbol', 'Strand': 'gene_strand',
                                        'Start': 'start_gene', 'End': 'end_gene',
                                        'Biotype': 'biotype', 'start': 'start_te',
                                        'end': 'end_te'})).sort_values(by=['chromosome', 'start_te']).
                        reset_index(drop=True))

    te_in_genes_sets = te_in_genes_sets[['chromosome', 'family', 'start_te', 'end_te', 'te_ins', 'no_te_ins',
                                         'gene_symbol', 'start_gene', 'end_gene', 'gene_strand', 'one_one', 'zero_zero',
                                         'log2FoldChange', 'lfcSE', 'biotype', 'start_bin', 'end_bin',
                                         'set_no', 'unique_no']]

    te_in_genes_sets['te_ins'] = clean_sample_list(te_in_genes_sets['te_ins'])
    te_in_genes_sets['no_te_ins'] = clean_sample_list(te_in_genes_sets['no_te_ins'])
    te_in_genes_sets['one_one'] = clean_sample_list(te_in_genes_sets['one_one'])
    te_in_genes_sets['zero_zero'] = clean_sample_list(te_in_genes_sets['zero_zero'])

    print(f'Sum of MITEs in DE genes and bins with the same sample sets: {len(te_in_genes_sets)}')
    print(te_in_genes_sets.to_string(max_rows=30))

    # te_in_genes_sets.to_csv('../files/P4_correct_MITEs_DE_bins_1.csv', sep='\t', index=False)
    # te_in_genes_sets.to_csv('../files/P1_MITEs_updown2000_DE_bins.csv', sep='\t', index=False)


# te_in_bins_and_genes(mite_list_p1, de_and_bins_p1)
# te_in_bins_and_genes(mite_list_p1, xloc_and_bins_p1)
# te_in_bins_and_genes(mite_list_p4, de_and_bins_p4)
te_in_bins_and_genes(mite_list_p4, xloc_and_bins_p4)

print("--- %s seconds ---" % (time.time() - start_time))
