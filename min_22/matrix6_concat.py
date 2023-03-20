import pandas as pd

df1 = pd.read_csv('family_6_matrix/hAT26_nonref_matrix.csv', sep='\t')
df1['family'] = 'hAT26'
df1ok = df1.reindex(columns=['chr', 'start', 'end', 'family', '1', '2', '3', '4', '5', '6', '7', '8', '9', '10', '11', '12'])

df2 = pd.read_csv('family_6_matrix/hAT27_nonref_matrix.csv', sep='\t')
df2['family'] = 'hAT27'
df2ok = df2.reindex(columns=['chr', 'start', 'end', 'family', '1', '2', '3', '4', '5', '6', '7', '8', '9', '10', '11', '12'])

df3 = pd.read_csv('family_6_matrix/hAT28_nonref_matrix.csv', sep='\t')
df3['family'] = 'hAT28'
df3ok = df3.reindex(columns=['chr', 'start', 'end', 'family', '1', '2', '3', '4', '5', '6', '7', '8', '9', '10', '11', '12'])

df4 = pd.read_csv('family_6_matrix/Stow3_nonref_matrix.csv', sep='\t')
df4['family'] = 'Stowaway3'
df4ok = df4.reindex(columns=['chr', 'start', 'end', 'family', '1', '2', '3', '4', '5', '6', '7', '8', '9', '10', '11', '12'])

df5 = pd.read_csv('family_6_matrix/Stow4_nonref_matrix.csv', sep='\t')
df5['family'] = 'Stowaway4'
df5ok = df5.reindex(columns=['chr', 'start', 'end', 'family', '1', '2', '3', '4', '5', '6', '7', '8', '9', '10', '11', '12'])

df6 = pd.read_csv('family_6_matrix/Stow36_nonref_matrix.csv', sep='\t')
df6['family'] = 'Stowaway36'
df6ok = df6.reindex(columns=['chr', 'start', 'end', 'family', '1', '2', '3', '4', '5', '6', '7', '8', '9', '10', '11', '12'])

df_final = pd.concat([df1ok, df2ok, df3ok, df4ok, df5ok, df6ok], ignore_index=True)
df_final_sort = df_final.sort_values(by=['chr', 'start'])

print(df_final_sort.head(), df_final_sort.tail())

df_final_sort.to_csv('matrix_nonref_6all.csv', sep='\t', index=False)
