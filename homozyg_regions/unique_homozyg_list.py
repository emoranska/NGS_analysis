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
print(df.to_string(max_rows=50), type(df.loc[0, 'genes_count']))

# create new df with unique sets of 1/1 and 0/0 samples
# df_unique = df.groupby(['1/1', '0/0']).size().reset_index().rename(columns={0: 'count'})
df_unique = df.groupby(['one_one_3first', 'zero_zero_3first'], as_index=False, sort=False).agg({'genes_count': 'sum'})

print(df_unique.to_string(max_rows=50))

# compare all bins with genes inside with unique list of samples sets and add the number of unique sets to bins
# (to be able to identify genes specific for every bins later)
'''
# method 1 - not ideal output
df['unique_no'] = None
for x in df.itertuples(index=True):
    for y in df_unique.itertuples(index=True):
        if x.one_one_3first == y.one_one_3first and x.zero_zero_3first == y.zero_zero_3first:
            df.at[x.Index, x.unique_no] = getattr(y, 'Index')

# method 2 (as 1)
for x in df.itertuples():
    for y in df_unique.itertuples():
        if x.one_one_3first == y.one_one_3first and x.zero_zero_3first == y.zero_zero_3first:
            df.at[x.Index, x.unique_no] = y.Index

# method 3 - very slow
for index_x, x in df.iterrows():
    for index_y, y in df_unique.iterrows():
        if x['one_one_3first'] == y['one_one_3first'] and x['zero_zero_3first'] == y['zero_zero_3first']:
            df.at[index_x, 'unique_no'] = index_y
'''
# method 4 - columns to merge on
cols = ["one_one_3first", "zero_zero_3first"]

# merge and take index, setting this as df["unique_no"]
df["unique_no"] = df[cols].merge(df_unique[cols].reset_index(), how="left")["index"]

print(df.to_string(max_rows=50))

# create column 'df_unique_index' to be able to keep row indexes after merging
df_unique = df_unique.reset_index(drop=False).rename(columns={'index': "df_unique_index"})

# create list of samples in 'one_zero' and 'zero_one' arrangement to be able to merge records with the same 1/1 and 0/0
# samples crosswise
df_unique['one_one_3first'] = df_unique['one_one_3first'].apply(ast.literal_eval)
df_unique['zero_zero_3first'] = df_unique['zero_zero_3first'].apply(ast.literal_eval)
df_unique['one_zero'] = (df_unique['one_one_3first'] + df_unique['zero_zero_3first']).apply(list)
df_unique['zero_one'] = (df_unique['zero_zero_3first'] + df_unique['one_one_3first']).apply(list)
df_unique['one_zero'] = df_unique['one_zero'].astype(str)
df_unique['zero_one'] = df_unique['zero_one'].astype(str)

print(df_unique.to_string(max_rows=50))

# merge records with the same 1/1 and 0/0 samples crosswise with keeping indexes of merged rows and counting the sum
# of genes
group = df_unique[['one_zero', 'zero_one']].apply(frozenset, axis=1)

df_unique_final_1 = df_unique.groupby(group).agg({'df_unique_index': (lambda x: list(x)), 'one_one_3first': 'first',
                                                  'zero_zero_3first': 'first','genes_count': 'sum'}).reset_index()

df_unique_final_1 = df_unique_final_1.drop(columns=['index']).rename(columns={'one_one_3first': 'one_one',
                                                                              'zero_zero_3first': 'zero_zero'}).\
    reindex(columns=['one_one', 'zero_zero', 'df_unique_index', 'genes_count'])

genes_sum = df_unique_final_1['genes_count'].sum()
print(df_unique_final_1.to_string(max_rows=50))
print(genes_sum)

# different solution tested
# df_unique_final3 = df_unique.groupby(group, as_index=False).first()
# df_unique_final = df_unique.groupby(group, as_index=False).agg({'genes_count': 'sum'})
# df_unique_final_1 = df_unique.groupby(group).sum('genes_count').reset_index()
# df_unique_final_1 = df_unique.groupby(group).agg({'one_one_3first': 'first', 'zero_zero_3first': 'first',
#                                                   'genes_count': 'sum'}).reset_index()
# genes_count_in_df_unique_final = df_unique.groupby(group, as_index=False, sort=False).agg({'genes_count': 'sum'})
# df_unique_final2 = df_unique.groupby(group, as_index=False).size().reset_index().rename(columns={0: 'count'})
# df_unique_final = df_unique_final.drop(columns=['genes_count', 'one_zero', 'zero_one'])
# df_unique_final['genes_count'] = genes_count_in_df_unique_final['genes_count']
# print(df_unique_final)
# print(df_unique_final.to_string(max_rows=50))
# print(df_unique_final3.to_string(max_rows=50))
# print(genes_count_in_df_unique_final)
# print(df_unique_final2.to_string(max_rows=50))

print("--- %s seconds ---" % (time.time() - start_time))
