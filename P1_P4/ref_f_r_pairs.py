import pandas as pd

df_raw = pd.read_csv('../P1_P4/NC_P1-6-to-95.teinsertions', sep="\t", decimal=',')

df_sort = df_raw.sort_values(by=['chr', 'pos', 'sample'])

just_f = df_sort[(df_sort["support"] == 'F')]
just_r = df_sort[(df_sort["support"] == 'R')]

# create merging df for F and R pairs and calculate mean for frequency F and R paired records --> to matrix
f_r_pairs = (just_f.merge(just_r, on=["sample", "chr", "strand", "family", "order", "comment"], suffixes=("_f", "_r")).
             query("pos_f <= pos_r <= pos_f + 1000")).sort_values(by=['chr', 'pos_f'])
f_r_pairs['freq_average'] = f_r_pairs[['frequency_f', 'frequency_r']].mean(axis=1).round(3)
print(f_r_pairs.to_string(max_rows=20))

f_from_pairs = f_r_pairs.drop(['pos_r', 'support_r', 'frequency_r'], axis=1)
f_from_pairs_to_concat = f_from_pairs.rename(columns={"pos_f":"pos", 'support_f':'support', 'frequency_f':"frequency"})
# print(f_from_pairs_to_concat.to_string(max_rows=20))

r_from_pairs = f_r_pairs.drop(['pos_f', "support_f", 'frequency_f'], axis=1)
r_from_pairs_to_concat = r_from_pairs.rename(columns={"pos_r":"pos", 'support_r':'support', 'frequency_r':"frequency"})
# print(r_from_pairs_to_concat.to_string(max_rows=20))

f_r_pairs_concat = pd.concat([f_from_pairs_to_concat, r_from_pairs_to_concat], ignore_index=True)
f_r_pairs_concat_sort = f_r_pairs_concat.sort_values(by=['chr', 'pos'])
f_r_pairs_concat_reset = f_r_pairs_concat_sort.reset_index(drop=True)
print(f_r_pairs_concat_reset.to_string(max_rows=20))

# saving file with F and R pairs one by one to prepare csv file with bins in Excel
f_r_pairs_concat_reset.to_csv("P1_F_R_pairs.csv", index=False, sep='\t')
