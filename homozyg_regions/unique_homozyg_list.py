import pandas as pd
import time
import ast

start_time = time.time()

df = pd.read_csv('bins/P1_EL10_chr1_9_homozyg_bins_with_genes.csv', sep='\t')

# create columns with 3 first samples for 1/1 and 0/0
df['1/1_3first'] = df['1/1'].apply(lambda x: ast.literal_eval(x)[:3])
df['0/0_3first'] = df['0/0'].apply(lambda x: ast.literal_eval(x)[:3])
df['genes_count'].astype(int)

df['1/1_3first'] = df['1/1_3first'].astype(str)
df['0/0_3first'] = df['0/0_3first'].astype(str)
# df['1/1_3first'] = df['1/1_3first'].apply(lambda x: ','.join(map(str, x)))
# df['0/0_3first'] = df['0/0_3first'].apply(lambda x: ','.join(map(str, x)))

df = df.rename(columns={'1/1_3first': 'one_one_3first', '0/0_3first': 'zero_zero_3first'}).drop(columns=['1/1', '0/0'])
df['unique_no'] = None
print(df.to_string(max_rows=50), type(df.loc[0, 'genes_count']))

# create new df with unique sets of 1/1 and 0/0 samples
# df_unique = df.groupby(['1/1', '0/0']).size().reset_index().rename(columns={0: 'count'})
df_unique = df.groupby(['one_one_3first', 'zero_zero_3first'], as_index=False, sort=False).agg({'genes_count': 'sum'})

print(df_unique.to_string(max_rows=50))

# compare all bins with genes inside with unique list of samples sets and add the number of unique sets to bins
# (to be able to identify genes specific for every bins later)

# for x in df.itertuples(index=True):
#     for y in df_unique.itertuples(index=True):
#         if x.one_one_3first == y.one_one_3first and x.zero_zero_3first == y.zero_zero_3first:
#             df.at[x.Index, x.unique_no] = getattr(y, 'Index')

# for x in df.itertuples():
#     for y in df_unique.itertuples():
#         if x.one_one_3first == y.one_one_3first and x.zero_zero_3first == y.zero_zero_3first:
#             df.at[x.Index, x.unique_no] = y.Index
#
# for index_x, x in df.iterrows():
#     for index_y, y in df_unique.iterrows():
#         if x['one_one_3first'] == y['one_one_3first'] and x['zero_zero_3first'] == y['zero_zero_3first']:
#             df.at[index_x, 'unique_no'] = index_y

print(df.to_string(max_rows=50))

df_unique['one_zero'] = df_unique['one_one_3first'] + df_unique['zero_zero_3first']
df_unique['zero_one'] = df_unique['zero_zero_3first'] + df_unique['one_one_3first']
print(df_unique.to_string(max_rows=50))

group = df_unique[['one_zero', 'zero_one']].apply(frozenset, axis=1)
df_unique_final = df_unique.groupby(group, as_index=False).first()
genes_count_in_df_unique_final = df_unique.groupby(group, as_index=False, sort=False).agg({'genes_count': 'sum'})
# df_unique_final2 = df_unique.groupby(group, as_index=False).size().reset_index().rename(columns={0: 'count'})
df_unique_final = df_unique_final.drop(columns=['genes_count', 'one_zero', 'zero_one'])
df_unique_final['genes_count'] = genes_count_in_df_unique_final['genes_count']
genes_sum = genes_count_in_df_unique_final['genes_count'].sum()

print(df_unique_final.to_string(max_rows=50))
# print(df_unique_final2.to_string(max_rows=50))
print(genes_count_in_df_unique_final)
print(genes_sum)
print("--- %s seconds ---" % (time.time() - start_time))
