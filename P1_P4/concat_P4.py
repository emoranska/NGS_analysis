import pandas as pd

b1 = pd.read_csv('../files/P4_MITEs_in_bins_1.csv', sep='\t')
b2 = pd.read_csv('../files/P4_MITEs_in_bins_2.csv', sep='\t')

b_all = pd.concat([b1, b2], ignore_index=True)
# b_all.to_csv('../files/P4_MITEs_in_bins_all.csv', sep='\t')

print(b_all.to_string(max_rows=30))

m1 = pd.read_csv('../files/P4_correct_MITEs_DE_bins_1.csv', sep='\t')
m2 = pd.read_csv('../files/P4_correct_MITEs_DE_bins_2.csv', sep='\t')

m_all = pd.concat([m1, m2], ignore_index=True)
# m_all.to_csv('../files/P4_correct_MITEs_DE_bins_all.csv', sep='\t')

unique_bins_count_m = m_all.drop_duplicates(subset=['start_bin', 'end_bin'])
print(f'Number of bins with MITEs in DE genes: {len(unique_bins_count_m)}')

print(m_all.to_string(max_rows=30))

ud1 = pd.read_csv('../files/P4_correct_MITEs_updown2000_DE_bins_1.csv', sep='\t')
ud2 = pd.read_csv('../files/P4_correct_MITEs_updown2000_DE_bins_2.csv', sep='\t')

ud_all = pd.concat([ud1, ud2], ignore_index=True)
ud_all.to_csv('../files/P4_correct_MITEs_updown2000_DE_bins_all.csv', sep='\t', index=False)

unique_bins_count_ud = ud_all.drop_duplicates(subset=['start_bin', 'end_bin'])
print(f'Number of bins with MITEs in DE genes: {len(unique_bins_count_ud)}')

print(ud_all.to_string(max_rows=30))

# add MITEs in XLOC genes to the MITEs list for annotated genes
xloc = pd.read_csv('../files/P4_correct_MITEs_XLOC_DE_bins.csv', sep='\t')

m_all_and_xloc = pd.concat([m_all, xloc], ignore_index=True).sort_values(by=['chromosome', 'start_te']).fillna('-')
print(m_all_and_xloc.to_string(max_rows=80))
m_all_and_xloc.to_csv('../files/P4_correct_MITEs_DE_bins_all_and_XLOC.csv', sep='\t', index=False)

xloc_ud = pd.read_csv('../files/P4_correct_MITEs_updown2000_XLOC_DE_bins.csv', sep='\t')
ud_all_and_xloc = (pd.concat([ud_all, xloc_ud], ignore_index=True).sort_values(by=['chromosome', 'start_te']).
                   fillna('-'))
print(ud_all_and_xloc.to_string(max_rows=50))
ud_all_and_xloc.to_csv('../files/P4_correct_MITEs_updown2000_DE_bins_all_and_XLOC.csv', sep='\t', index=False)

# check if there are the same MITEs in P1 and P4
m_p1 = pd.read_csv('../files/P1_correct_MITEs_DE_bins.csv', sep='\t')
check_the_same_mites = m_all.merge(m_p1, how='inner', on=['chromosome', 'family', 'start_te', 'end_te'])
print(check_the_same_mites.to_string())
