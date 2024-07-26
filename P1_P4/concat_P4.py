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
