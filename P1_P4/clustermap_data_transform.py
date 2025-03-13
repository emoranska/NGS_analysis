import pandas as pd
import ast
import numpy as np

df_matrix = pd.read_csv('../files/clustermaps/P4_DE005_cm_general_data.csv', sep='\t')
print(df_matrix.to_string())

df_matrix['one_one'] = df_matrix['one_one'].apply(lambda x: ast.literal_eval(x))
df_matrix['zero_zero'] = df_matrix['zero_zero'].apply(lambda x: ast.literal_eval(x))
print(df_matrix.to_string())

df_matrix_trans = df_matrix.assign(
    **{f"{k}_{i+1}": df_matrix.values[
    np.arange(len(df_matrix)),
    df_matrix.columns.get_indexer(df_matrix[k].str[i])]
       for k in ['one_one', 'zero_zero'] for i in range(3)})

print(df_matrix_trans.to_string(max_rows=30))