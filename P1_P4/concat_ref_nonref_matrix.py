import pandas as pd

ref = pd.read_csv('../min_22/1-DH-to-12-RO_ref_matrix_sort_2000_all.csv', sep="\t")
nonref = pd.read_csv('../min_22/1-DH-to-12-RO_nonref_matrix_all.csv', sep="\t")

# ref = pd.read_csv('../P1_P4/P1_ref_matrix_sort_2000_all.csv', sep="\t")
# nonref = pd.read_csv('../P1_P4/P1_nonref_matrix_all.csv', sep="\t")

matrix_final = pd.concat([ref, nonref], ignore_index=True).sort_values(by=['chr', 'start'])

matrix_final.to_csv('../min_22/1-DH-to-12-RO_matrix_all.csv', index=False, sep='\t')
print(matrix_final.to_string(max_rows=20))
