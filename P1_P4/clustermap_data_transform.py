import pandas as pd
import ast
import numpy as np

df_matrix = pd.read_csv('../files/heatmaps/cm_data/P4_all_genes_with_MITEs_cm_general_data.csv', sep='\t')

print(df_matrix.to_string())

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
print(df_matrix_trans.to_string(max_rows=30))

df_final = (df_matrix_trans.rename(columns={'one_one_1': 'ins_1', 'one_one_2': 'ins_2', 'one_one_3': 'ins_3',
                                           'zero_zero_1': 'no_ins_1', 'zero_zero_2': 'no_ins_2',
                                           'zero_zero_3': 'no_ins_3'}).
            drop(columns=df_matrix_trans.loc[:, 'one_one':'no_te_ins']).
            drop(columns=df_matrix_trans.loc[:, 'ins_up_1': 'ins_down'])).drop(columns=['de_results'])
print(df_final.to_string(max_rows=30))
df_final.to_csv('../files/heatmaps/cm_data/P4_all_genes_with_MITEs_cm_general_data_final.csv', sep='\t', index=False)
