import pandas as pd

te_final = pd.read_csv('../files/P4_MITEs_DE_bins.csv', sep='\t')
te_ref = pd.read_csv('../files/P4_ref_matrix_sort_2000_all.csv', sep='\t')
te_nonref = pd.read_csv('../files/P4_nonref_matrix_all.csv', sep='\t')

ref_check = te_final.merge(te_ref, on=['chr', 'family', 'start_te', 'end_te'])
print(f'Reference MITES: {ref_check}')

nonref_check = te_final.merge(te_nonref, on=['chr', 'family', 'start_te', 'end_te'])
print(f'Non-reference MITES: {nonref_check}')
