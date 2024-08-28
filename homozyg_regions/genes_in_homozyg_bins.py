import pandas as pd
import time
import ast

start_time = time.time()

# df = pd.read_csv('~/Pulpit/burak/homozyg_regions/bins/P4_EL10_chr1_9_homozyg_bins_with_genes.csv', sep='\t')
df = pd.read_csv('~/Pobrane/P1_rest_bins_with_genes.csv', sep='\t')

# create columns with 3 first samples for 1/1 and 0/0
df['1/1_3first'] = df['1/1'].apply(lambda x: ast.literal_eval(x)[:3])
df['0/0_3first'] = df['0/0'].apply(lambda x: ast.literal_eval(x)[:3])
df['genes_count'].astype(int)

df['1/1_3first'] = df['1/1_3first'].astype(str)
df['0/0_3first'] = df['0/0_3first'].astype(str)

df = df.rename(columns={'1/1_3first': 'one_one_3first', '0/0_3first': 'zero_zero_3first'}).drop(columns=['1/1', '0/0'])

# create new df with unique sets of 1/1 and 0/0 samples
df_unique = df.groupby(['one_one_3first', 'zero_zero_3first'], as_index=False, sort=False).agg({'genes_count': 'sum'})

print(df_unique.to_string(max_rows=50))

# compare all bins with genes inside with unique list of samples sets and add the number of unique sets to bins
# (to be able to identify genes specific for every bins later)
# method 4 - columns to merge on
cols = ["one_one_3first", "zero_zero_3first"]

# merge and take index, setting this as df["unique_no"]
# df["unique_no"] = df[cols].merge(df_unique[cols].reset_index(), how="left")["index"]
print(df.to_string(max_rows=50))

# df_with_genes = pd.read_csv('~/Pulpit/burak/homozyg_regions/bins/P4_chr1_9_genes_in_homozyg_bins.csv',
#                             sep='\t')
df_with_genes = pd.read_csv('~/Pobrane/P1_chr1_9_genes_in_rest_bins_to_add.csv', sep='\t')
print(df_with_genes.to_string(max_rows=50))

df_merged = df_with_genes.merge(df, on=['chr', 'start', 'end']).rename(columns={'start': 'start_bin', 'end': 'end_bin',
                                                                                'one_one_3first': 'one_one',
                                                                                'zero_zero_3first': 'zero_zero'}).\
    drop(columns=['chr_gene']).reindex(columns=['chr', 'start_bin', 'end_bin', 'start_gene', 'end_gene', 'gene_name',
                                                'gene_strand', 'unique_no', 'one_one', 'zero_zero', 'genes_count'])
print(df_merged.to_string(max_rows=50))

# df_merged.to_csv('P4_list_genes_in_homozyg_bins.csv', sep='\t', index=False)
df_merged.to_csv('../P1_list_genes_in_rest_bins_to_add.csv', sep='\t', index=False)
