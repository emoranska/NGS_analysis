import pandas as pd
import numpy as np

df_raw = pd.read_csv('1-DH-to-12-RO.NC_teinsertions.txt', sep="\t", decimal=',')
# print(df.head(5))

df_sort = df_raw.sort_values(by=['chr', 'pos', 'sample'])

family = "Stowaway|4"

non_ref = df_sort[(df_sort["family"] == family) & (df_sort["support"] == 'FR')]
non_ref_reset = non_ref.reset_index(drop=True)
ref = df_sort[(df_sort["family"] == family) & (df_sort["support"] != 'FR')]

# non_ref.to_csv("Stowaway4_FR.csv", index=False, sep='\t')
# ref.to_csv("Stowaway4_F_R.csv", index=False, sep='\t')

just_f = df_sort[(df_sort["family"] == family) & (df_sort["support"] == 'F')]
just_r = df_sort[(df_sort["family"] == family) & (df_sort["support"] == 'R')]
# print(just_f.head())
# print(just_r.head())

f_r_pairs = pd.DataFrame(columns=just_f.columns)

# choosing rows for reference TE insertions (having pairs with F and R in range 1000 bp)
for f in just_f.itertuples():
    for r in just_r.itertuples():
        if f.sample == r.sample and f.chr == r.chr and r.pos in range(f.pos, f.pos + 1000):
            f_r_pairs = f_r_pairs.append(pd.DataFrame([f]))
            f_r_pairs = f_r_pairs.append(pd.DataFrame([r]))
# print(f_r_pairs)

f_r_pairs_drop = f_r_pairs.drop(['Index'], axis=1)

# choose unpaired records with F or R
unpaired_f_r = ref[~ref.apply(tuple, 1).isin(f_r_pairs_drop.apply(tuple, 1))]
unpaired_f_r_reset = unpaired_f_r.reset_index(drop=True)
# print(unpaired_f_r_reset.to_string())

'''
# preparing file for creating bins for reference insertions with bedtools merge
f_r_pairs_sort = f_r_pairs.sort_values(by=['chr', 'pos'])
f_r_pairs_bins = f_r_pairs_sort.drop(["strand", "family", "order", "support", "comment", "frequency", "Index"], axis=1)
f_r_pairs_bins['end'] = f_r_pairs_bins.loc[:, 'pos']
f_r_pairs_bins_final = f_r_pairs_bins.reindex(columns=['chr', 'pos', 'end', 'sample'])
f_r_pairs_bins_final.to_csv('Stow36_ref_bins.bed',  sep='\t', index=False, header=False)
'''

matrix_ref = pd.read_csv('Stow4_ref_merged1000.bed', sep='\t')

f_r_pairs_reset_drop = f_r_pairs_drop.reset_index(drop=True)

# another solutions for finding mean of paired reference rows (with F and R)

# ref_mean = f_r_pairs_reset_drop.groupby(f_r_pairs_reset_drop.index).mean()
# f_r_pairs_reset_drop["frequency"] = f_r_pairs_reset_drop["frequency"].astype(float)
# f_r_pairs_reset_drop["freq_average"] = f_r_pairs_reset_drop["frequency"].rolling(2).mean()

# f_r_pairs_reset_drop.frequency = pd.to_numeric(f_r_pairs_reset_drop.frequency, errors='coerce')
# f_r_pairs_mean = f_r_pairs_reset_drop.groupby(np.arange(f_r_pairs_reset_drop.shape[0]) // 2).agg({'sample': 'sample', 'frequency': 'mean'}).fillna(0)

# merge paired TE reference records (with F and R for the sample) and calculate mean from frequency of F and R
f_r_pairs_mean = f_r_pairs_reset_drop.assign(pos_end=f_r_pairs_reset_drop.pos.shift(-1), freq_2=f_r_pairs_reset_drop.frequency.shift(-1))[::2]
f_r_pairs_mean['freq_average'] = f_r_pairs_mean[['frequency', 'freq_2']].mean(axis=1).round(3)
# print(f_r_pairs_mean.to_string())

# creating matrix for reference insertions
for ref in f_r_pairs_mean.itertuples():
    for x in matrix_ref.itertuples():
        if ref.chr == x.chr and ref.pos in range(int(x.start), int(x.end)) and ref.pos_end in range(int(x.start), int(x.end)):
            matrix_ref.at[x.Index, ref.sample] = ref.freq_average

# print(matrix_ref.to_string())

matrix_ref.columns = matrix_ref.columns.astype("str")
matrix_ref_drop = matrix_ref.drop(['no_of_samples'], axis=1)
matrix_ref_final = matrix_ref_drop.reindex(columns=['chr', 'start', 'end', '1', '2', '3', '4', '5', '6', '7', '8', '9', '10', '11', '12'])
matrix_ref_final_2 = matrix_ref_final.fillna(0)
matrix_ref_final_2.to_csv('Stow4_ref_matrix.csv', sep='\t')
# print(matrix_ref_final_2.to_string())

nonref_unpaired = pd.DataFrame(columns=non_ref_reset.columns)
# print(nonref_unpaired)

for fr in non_ref_reset.itertuples():
    for un in unpaired_f_r_reset.itertuples():
        if un.sample != fr.sample and un.chr == fr.chr and un.pos in range(fr.pos - 200, fr.pos + 200):
            nonref_unpaired = nonref_unpaired.append(pd.DataFrame([un]))
# print(nonref_unpaired.to_string())

nonref_unpaired_final = nonref_unpaired.drop_duplicates()
# print(nonref_unpaired_final.to_string())

nonref_all = pd.concat([non_ref_reset, nonref_unpaired_final], ignore_index=True)
nonref_all_drop = nonref_all.drop(['Index'], axis=1)
nonref_all_drop_reset = nonref_all_drop.reset_index(drop=True)
nonref_all_sort = nonref_all_drop_reset.sort_values(by=['chr', 'pos'])
# print(nonref_all_sort.to_string())

matrix_nonref = pd.read_csv('Stow4_FRmerged200.bed', sep='\t')

# creating matrix for nonref TE insertions
for te in nonref_all_sort.itertuples():
    for x in matrix_nonref.itertuples():
        if te.chr == x.chr and te.pos in range(int(x.start), int(x.end+1)):
            matrix_nonref.at[x.Index, te.sample] = te.frequency

'''
# creating matrix for nonref TE insertions --> working for hAT, not for Stowaway (KeyError: 12)
for te in non_ref_reset.itertuples():
    for te_un in unpaired_f_r_reset.itertuples():
        for x in matrix_nonref.itertuples():
            if te.chr == x.chr and te.pos in range(int(x.start), int(x.end+1)):  # and te.sample in cols:
                matrix_nonref.at[x.Index, te.sample] = te.frequency
            elif te_un.chr == x.chr and te_un.pos in range(int(x.start - 100), int(x.end + 100)) and float(matrix_nonref.at[x.Index, te.sample]) < te_un.frequency:
                matrix_nonref.at[x.Index, te_un.sample] = te_un.frequency
'''

matrix_nonref.columns = matrix_nonref.columns.astype("str")
# print(matrix_nonref.at[2, '1'], type(matrix_nonref.at[2, '1']), pd.isnull(matrix_nonref.at[2, '1']))

matrix_nonref_drop = matrix_nonref.drop(['no_of_samples'], axis=1)
matrix_nonref_final = matrix_nonref_drop.reindex(columns=['chr', 'start', 'end', '1', '2', '3', '4', '5', '6', '7', '8', '9', '10', '11', '12'])
matrix_nonref_final_2 = matrix_nonref_final.fillna(0)
# print(matrix_nonref_final_2.to_string())
matrix_nonref_final_2.to_csv('Stow4_nonref_matrix.csv', sep='\t')
