import pandas as pd

genes_bins = pd.read_csv('../files/P1_list_genes_in_homozyg_bins.csv', sep='\t')
genes_bins_rest = pd.read_csv('../files/P1_list_genes_in_rest_bins_to_add.csv', sep='\t')

genes_bins['bin_length'] = genes_bins['end_bin'] - genes_bins['start_bin']
col = genes_bins.pop('bin_length')
genes_bins.insert(3, col.name, col)
print(genes_bins.to_string(max_rows=30))

genes_bins_rest['bin_length'] = genes_bins_rest['end_bin'] - genes_bins_rest['start_bin']
col = genes_bins_rest.pop('bin_length')
genes_bins_rest.insert(3, col.name, col)
print(genes_bins_rest.to_string(max_rows=30))

all_genes_bins = ((pd.concat([genes_bins, genes_bins_rest], ignore_index=True).sort_values(by=['chr', 'start_bin']).
                  drop_duplicates(subset=['chr', 'start_gene', 'end_gene', 'gene_name', 'unique_no'], keep='last')).
                  reset_index(drop=True))
print(all_genes_bins.to_string(max_rows=500))

genes_count = all_genes_bins.drop_duplicates(subset=['gene_name'])
print(len(genes_count))

sets_list = pd.read_csv('../files/P1_list_to_DE_20more_genes.csv', sep='\t')
all_genes_bins = (sets_list.merge(all_genes_bins, on=['unique_no']).
                  drop(columns=['genes_count_x', 'df_unique_index', 'genes_count_y']).reset_index(drop=True))

# set_no = all_genes_bins.pop('set_no')
# all_genes_bins.insert(8, set_no.name, set_no)

# all_genes_bins.to_csv('../files/P1_all_EL10_genes_in_all_bins.csv', sep='\t', index=False)

print(all_genes_bins.to_string(max_rows=50))
