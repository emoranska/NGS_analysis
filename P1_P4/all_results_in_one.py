import pandas as pd

bins_and_de = pd.read_csv('../files/P1_all_EL10_genes_in_all_bins_with_DE_and_kallisto.csv', sep='\t',
                          dtype={'log2FoldChange': 'float64', 'lfcSE': 'float64'})
mites = pd.read_csv('../files/P1_correct_MITES_in_genes_and_updown.csv', sep='\t',
                    dtype={'log2FoldChange': 'float64', 'lfcSE': 'float64'})
print(mites.to_string(max_rows=30))

bins_and_de = (bins_and_de.drop(columns=['start_gene', 'end_gene']).
               rename(columns={'Start': 'start_gene', 'End': 'end_gene', 'gene_name': 'gene_symbol'}))

start_gene = bins_and_de.pop('start_gene')
bins_and_de.insert(5, start_gene.name, start_gene)
bins_and_de['start_gene'] = bins_and_de['start_gene']

end_gene = bins_and_de.pop('end_gene')
bins_and_de.insert(6, end_gene.name, end_gene)
bins_and_de['end_gene'] = bins_and_de['end_gene']

print(bins_and_de.to_string(max_rows=50))

all_together = (((bins_and_de.merge(mites, how='outer', on=['chr', 'start_bin', 'end_bin', 'bin_length',
                                                         'gene_symbol', 'start_gene', 'end_gene', 'gene_strand',
                                                         'set_no', 'unique_no']).
                drop(columns=['biotype', 'gene_function', 'one_one_y', 'zero_zero_y'])).
                rename(columns={'one_one_x': 'one_one', 'zero_zero_x': 'zero_zero'})).
                sort_values(by=['chr', 'start_gene'])).reset_index(drop=True)

print(all_together.to_string(max_rows=200))
