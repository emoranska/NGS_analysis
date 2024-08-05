import pandas as pd

te_final = pd.read_csv('../files/P1_MITEs_updown2000_DE_bins.csv', sep='\t')
te_ref = pd.read_csv('../files/P1_ref_matrix_sort_2000_all.csv', sep='\t')
te_nonref = pd.read_csv('../files/P1_nonref_matrix_all.csv', sep='\t')

ref_check = te_final.merge(te_ref, on=['chr', 'family', 'start_te', 'end_te'])
print(f'Reference MITEs: \n {ref_check}')

nonref_check = te_final.merge(te_nonref, on=['chr', 'family', 'start_te', 'end_te'])
print(f'Non-reference MITEs: \n {nonref_check}')
