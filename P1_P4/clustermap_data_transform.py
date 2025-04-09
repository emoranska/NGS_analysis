import pandas as pd
import ast
import numpy as np

df_matrix = pd.read_csv('../files/heatmaps/cm_data/P1_DE005_with_MITEs_cm_general_data.csv', sep='\t')
# df_matrix = pd.read_csv('../files/heatmaps/cm_data/P1_all_genes_with_MITEs_cm_general_data_1.csv', sep='\t')

print(df_matrix.to_string(max_rows=30))

df_matrix['one_one'] = df_matrix['one_one'].apply(lambda x: ast.literal_eval(x))
df_matrix['zero_zero'] = df_matrix['zero_zero'].apply(lambda x: ast.literal_eval(x))
print(df_matrix.to_string(max_rows=30))

df_matrix_trans = df_matrix.assign(**{f"{k}_{i+1}": df_matrix.values[np.arange(
    len(df_matrix)), df_matrix.columns.get_indexer(df_matrix[k].str[i])]
       for k in ['one_one', 'zero_zero'] for i in range(3)})

print(df_matrix_trans.to_string(max_rows=30))

m = df_matrix_trans['ins_de'] != df_matrix_trans['de_results']
df_matrix_trans.loc[m, ['one_one_1', 'zero_zero_1']] = df_matrix_trans.loc[
    m, ['zero_zero_1', 'one_one_1']].values  # Swap rows if condition met

df_matrix_trans.loc[m, ['one_one_2', 'zero_zero_2']] = df_matrix_trans.loc[
    m, ['zero_zero_2', 'one_one_2']].values

df_matrix_trans.loc[m, ['one_one_3', 'zero_zero_3']] = df_matrix_trans.loc[
    m, ['zero_zero_3', 'one_one_3']].values

print(type(df_matrix_trans), '\n', df_matrix_trans.to_string(max_rows=30))

df_final = df_matrix_trans.rename(columns={'one_one_1': 'ins_1', 'one_one_2': 'ins_2', 'one_one_3': 'ins_3',
                                           'zero_zero_1': 'no_ins_1', 'zero_zero_2': 'no_ins_2',
                                           'zero_zero_3': 'no_ins_3'})
#            drop(columns=df_matrix_trans.loc[:, 'one_one':'no_te_ins']))
#             drop(columns=df_matrix_trans.loc[:, 'ins_up_1': 'ins_down'])).drop(columns=['de_results', 'ins_check'])
df_final = df_final.loc[:, ['gene_name', 'ins_1', 'ins_2', 'ins_3', 'no_ins_1', 'no_ins_2', 'no_ins_3', 'mite_loc',
                            'ins_de', 'family']]
print(df_final.to_string(max_rows=30))

df_final = (df_final.drop_duplicates().drop_duplicates(subset=['gene_name', 'mite_loc', 'ins_de', 'family']).
            reset_index(drop=True))
print(df_final.to_string(max_rows=30))

# Series as name for shorter reference
s = df_final['gene_name']
# group consecutive occurrences
group = s.ne(s.shift()).cumsum()
# form group and save as "g" for efficiency
g = s.groupby(group)
# identify groups with more than 1 value
m = g.transform('size').gt(1)
# increment values
df_final.loc[m, 'gene_name'] += '_TE'+g.cumcount().add(1).astype(str)
print(df_final.to_string(max_rows=30))
df_final.to_csv('../files/heatmaps/cm_data/P1_DE005_with_MITEs_cm_general_data_final_1.csv', sep='\t', index=False)
# df_final.to_csv('../files/heatmaps/cm_data/P1_all_genes_with_MITEs_cm_general_data_final_1.csv', sep='\t', index=False)
