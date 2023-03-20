import pandas as pd

df_raw = pd.read_csv('NC_gDNA-P4-1-to-64.teinsertions', sep="\t", decimal=',')

df_sort = df_raw.sort_values(by=['chr', 'pos', 'sample'])
'''
# FIRST PART code to run, then comment
non_ref = df_sort[(df_sort["support"] == 'FR')]
non_ref_reset = non_ref.reset_index(drop=True)
print(non_ref_reset.to_string(max_rows=20))
non_ref_reset.to_csv('../P1_P4/P4_FR.csv', index=False, sep='\t')
'''

# SECOND PART code to run
ref = df_sort[(df_sort["support"] != 'FR')]

f_r_pairs_reset_drop = pd.read_csv('P4_F_R_pairs.csv', sep="\t")

# choose unpaired records with F or R
unpaired_f_r = ref[~ref.apply(tuple, 1).isin(f_r_pairs_reset_drop.apply(tuple, 1))]
unpaired_f_r_reset = unpaired_f_r.reset_index(drop=True)
# print(unpaired_f_r_reset.to_string(max_rows=20))
# unpaired_f_r_reset.to_csv('P1_unpaired_F_R.csv', index=False, sep='\t')

non_ref_reset = pd.read_csv('P4_FR_chr/P4_FR_chr9.csv', sep='\t')
non_ref_reset_tomerge = non_ref_reset.drop(['strand', 'order', 'support', 'comment', 'frequency'], axis=1)
print(non_ref_reset_tomerge.to_string(max_rows=20))

# unpaired_f_r_reset = pd.read_csv('unpaired_F_R.csv', sep='\t')

# choose records from unpaired F and R which will be joined as a supplement to FR nonref records
nonref_unpaired = (non_ref_reset_tomerge.merge(unpaired_f_r_reset, on=["chr", "family"], suffixes=("_fr", "_un")).
                   query("pos_fr - 200 <= pos_un <= pos_fr + 200")).drop('pos_fr', 1).\
    query("sample_fr != sample_un").drop('sample_fr', 1).sort_values(by=['chr', 'pos_un']).\
    drop_duplicates().rename(columns={'sample_un':'sample', 'pos_un':'pos'}).\
    reindex(columns=['sample', 'chr', 'pos', 'strand', 'family', 'order', 'support', 'comment', 'frequency'])

print(nonref_unpaired.to_string(max_rows=20))
nonref_unpaired.to_csv('../P1_P4/P4_nonref_unpaired_chr/P4_nonref_unpaired_chr9.csv', sep='\t', index=False)
