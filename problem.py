import pandas as pd
import numpy as np
import time
import matplotlib.pyplot as plt
import seaborn as sns

start_time = time.time()

# d1 = {'chr': [1, 1], 'start': [64, 1000], 'end': [150, 2000], 'family': ['a', 'b'],
#       'ins': [['P1-12', 'P1-22', 'P1-25', 'P1-28', 'P1-90'],
#               ['P1-6', 'P1-89', 'P1-92', 'P1-93']],
#       'no_ins': [['P1-6', 'P1-89', 'P1-92', 'P1-93'],
#                  ['P1-12', 'P1-25', 'P1-28', 'P1-89', 'P1-90', 'P1-93']]}
# df1 = pd.DataFrame.from_dict(data=d1)
# print(df1.to_string())
#
# d2 = {'set_1': [['P1-12', 'P1-25', 'P1-28'], ['P1-6', 'P1-89', 'P1-93']],
#       'set_2': [['P1-89', 'P1-92', 'P1-93'], ['P1-25', 'P1-28', 'P1-90']]}
# df2 = pd.DataFrame.from_dict(data=d2)
# print(df2.to_string())
#
# d3 = {'chr': [1, 1], 'start': [64, 64], 'end': [150, 150], 'family': ['a', 'a'], 'df2_index': [0, 0],
#       'ins_set1': [['P1-12', 'P1-25', 'P1-28',], []],
#       'ins_set2': [[], ['P1-25', 'P1-89', 'P1-90']],
#       'no_ins_set1': [[], ['P1-12', 'P1-6', 'P1-92']],
#       'no_ins_set2': [['P1-89', 'P1-92', 'P1-93'], []]
#       }
# df3 = pd.DataFrame.from_dict(data=d3)
# print(df3.to_string())

# d1 = {'chr': [1, 1], 'start': [64, 1000], 'end': [150, 2000], 'family': ['a', 'b'],
#       'ins': [['P1-12', 'P1-22', 'P1-25', 'P1-28', 'P1-90'],
#               ['P1-6', 'P1-89', 'P1-92', 'P1-93']]}
# df1 = pd.DataFrame.from_dict(data=d1)
# print(df1.to_string())
#
# d2 = {'set_1': [['P1-12', 'P1-25', 'P1-28'], ['P1-6', 'P1-89', 'P1-93']],
#       'set_2': [['P1-89', 'P1-92', 'P1-93'], ['P1-25', 'P1-28', 'P1-90']]}
# df2 = pd.DataFrame.from_dict(data=d2)
# print(df2.to_string())
#
# d3 = {'chr': [1, 1, 1, 1], 'start': [64, 64, 1000, 1000], 'end': [150, 150, 2000, 2000], 'family': ['a', 'a', 'b', 'b'],
#       'df2_index': [0, 1, 0, 1],
#       'ins_set1': [['P1-12', 'P1-25', 'P1-28',], [], [], ['P1-6', 'P1-89', 'P1-93']],
#       'ins_set2': [[], ['P1-25', 'P1-28', 'P1-90'], ['P1-89', 'P1-92', 'P1-93'], []]}
# df3 = pd.DataFrame.from_dict(data=d3)
# print(df3.to_string())
#
# matches = [x for x in df2.iloc[0, 0] if x in df1.iloc[0, 4]]
# print(matches)
#
# d_sets = {'gene_id': ['g0', 'g2', 'g3'],
#           'set_1': [['P1-12', 'P1-25', 'P1-28',], ['P1-6', 'P1-12', 'P1-25'], ['P1-6', 'P1-25', 'P1-28']]}
#
# df_sets = pd.DataFrame.from_dict(data=d_sets)
# print(df_sets.to_string())
#
# d_counts = {'gene_id': ['g0', 'g1', 'g2', 'g3'],
#             'P1-6': [1, 2, 5, 7],
#             'P1-12': [2, 4, 6, 5],
#             'P1-25': [5, 1, 2, 6],
#             'P1-28': [3, 8, 1, 2],
#             }
# df_counts = pd.DataFrame.from_dict(data=d_counts)
# print(df_counts.to_string())
#
# # df_sets[['1_1', '1_2', '1_3']] = pd.DataFrame(df_sets.set_1.tolist(), index=df_sets.index)
#
# print(df_sets.to_string())
#
# # df = df_sets.merge(df_counts, on='gene_id')
# # print(df.to_string())
#
# out = df_sets.merge(df_counts, on='gene_id', how='left')
# print(out)
# m = out['set_1'].str.join('|').str.get_dummies().astype(bool)
# print(m)
# out.loc[:, m.columns] = out.loc[:, m.columns].where(m)
#
# out['mean_set_1'] = out.loc[:, m.columns].mean(axis=1)
#
# print(out)
#
#
# rpk_data = {'P1-6': [1, 2, 3, 4],
#             'P1-12': [2, 4, 6, 8],
#             'P1-25': [6, 12, 3, 9]
#             }
# rpk = pd.DataFrame.from_dict(rpk_data)
# print(rpk.to_string())
#
# sc_f_data = {'sample': ['P1-6', 'P1-12', 'P1-25'],
#              'factor': [1, 2, 3]
#              }
# scaling_factor = pd.DataFrame.from_dict(sc_f_data)
# print(scaling_factor.to_string())
#
# divided_value = rpk.loc[0, 'P1-6'] / scaling_factor.loc[0, 'factor']
# print(divided_value)
#
# out = rpk / scaling_factor.set_index('sample')['factor']
# print(out.to_string())
#
#
# plot_data = {'P1-6': [1, 2, 3], 'P1-12': [6, 8, 10], 'P1-25': [5, 2, 7]}
# df_plot = pd.DataFrame.from_dict(plot_data)
#
# # działa
# # plt.imshow(df_plot, cmap="RdYlBu")
#
# # działa
# fig, ax = plt.subplots(figsize=(12, 7))
# sns.heatmap(df_plot, cmap="YlGnBu", annot=True, fmt="0.0f")
# plt.show()

# # działa
# plt.pcolor(df_plot)
# # plt.yticks(np.arange(0.5, len(rpk_data.index), 1), rpk_data.index)
# # plt.xticks(np.arange(0.5, len(rpk_data.columns), 1), rpk_data.columns)
# plt.show()

df_matrix_data = {'11': [['P4-1', 'P4-2', 'P4-3'], ['P4-1', 'P4-3', 'P4-4']],
                  '00': [['P4-4', 'P4-6', 'P4-7',], ['P4-2', 'P4-5', 'P4-7']],
                  'P4-1': [1, 2], 'P4-2': [6, 8], 'P4-3': [5, 2], 'P4-4': [2, 3], 'P4-5': [np.nan, 2], 'P4-6': [6, np.nan],
                  'P4-7': [3, 2]}
df_matrix = pd.DataFrame.from_dict(df_matrix_data)
print(df_matrix.to_string())

df_output_data = {'11': [['P4-1', 'P4-2', 'P4-3'], ['P4-1', 'P4-3', 'P4-4']],
                  '00': [['P4-4', 'P4-6', 'P4-7',], ['P4-2', 'P4-5', 'P4-7']],
                  'P4-1': [1, 2], 'P4-2': [6, 8], 'P4-3': [5, 2], 'P4-4': [2, 3], 'P4-5': [np.nan, 2], 'P4-6': [6, np.nan],
                  'P4-7': [3, 2], '11_1': [1, 2], '11_2': [6, 2], '11_3': [5, 3],
             '00_1': [2, 8], '00_2': [6, 2], '00_3': [3, 2]}

df_output = pd.DataFrame.from_dict(df_output_data)
print(df_output.to_string())

# def extract_values(row, col_name):
#     return [row[item] if item in row else np.nan for item in row[col_name]]
#
# for col in ['11', '00']:
#     extracted_values = df_matrix.apply(lambda row: extract_values(row, col), axis=1)
#     df_expanded = pd.DataFrame(extracted_values.tolist(), columns=[f"{col}_{i+1}" for i in range(extracted_values.str.len().max())])
#     df_matrix = pd.concat([df_matrix, df_expanded], axis=1)

# # columns to consider
# cols = ['11', '00']

# # first reshape the columns with the lists
# tmp1 = (df_matrix[cols]
#         .melt(value_name='col', ignore_index=False)
#         .explode('col').reset_index()
#         .assign(n=lambda x: x.groupby(['index', 'variable']).cumcount()+1)
#        )
# # then reshape the columns with the values
# tmp2 = (df_matrix.drop(columns=cols, errors='ignore')
#         .melt(var_name='col', ignore_index=False)
#         .reset_index()
#        )
#
# # merge, reshape, rename columns
# out = tmp1.merge(tmp2, how='left').pivot(index='index', columns=['variable', 'n'], values='value')
# out.columns = out.columns.map(lambda x: f'{x[0]}_{x[1]}')
#
# # join to original
# out = df_matrix.join(out)
# print(out.to_string())

# df_matrix_trans = df_matrix.assign(
#     **{f"{k}_{i+1}": df_matrix.values[
#     np.arange(len(df_matrix)),
#     df_matrix.columns.get_indexer(df_matrix[k].str[i])]
#        for k in ['11', '00'] for i in range(3)})
#
# print(df_matrix_trans.to_string())
# print("--- %s seconds ---" % (time.time() - start_time))

df_genes_data = {'gene_id': ['g0', 'g1', 'g1', 'g2', 'g3', 'g4', 'g4', 'g4']}
df_genes = pd.DataFrame.from_dict(df_genes_data)
print(df_genes.to_string())

df_genes_marked_data = {'gene_id': ['g0', 'g1_TE1', 'g1_TE2', 'g2', 'g3', 'g4_TE1', 'g4_TE2', 'g4_TE3']}
df_genes_marked = pd.DataFrame.from_dict(df_genes_marked_data)
print(df_genes_marked.to_string())

rep = []
gene_list = df_genes['gene_id']
for idx in range(0, len(gene_list) - 1):
    # getting Consecutive elements
    if gene_list[idx] == gene_list[idx + 1]:
        rep.append(gene_list[idx])

# getting list of repetitive gene names
rep = list(set(rep))

# printing result
print("Consecutive identical elements are : " + str(rep))

# Series as name for shorter reference
s = df_genes['gene_id']
# group consecutive occurrences
group = s.ne(s.shift()).cumsum()
# form group and save as "g" for efficiency
g = s.groupby(group)
# identify groups with more than 1 value
m = g.transform('size').gt(1)
# increment values
df_genes.loc[m, 'gene_id'] += '_TE'+g.cumcount().add(1).astype(str)

print(df_genes.to_string())