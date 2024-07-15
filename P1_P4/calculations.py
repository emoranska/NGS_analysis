import pandas as pd
import time

start_time = time.time()

te_matrix = pd.read_csv('../files/P1_matrix_all.csv', sep='\t')
# te_matrix = pd.read_csv('../files/P4_matrix_all.csv', sep='\t')
print(te_matrix.to_string(max_rows=30))

te_matrix['start'] = te_matrix['start'].astype(int)
te_matrix['end'] = te_matrix['end'].astype(int)

cols = te_matrix.filter(like='P1-').columns
# cols = te_matrix.filter(like='P4-').columns
sep = ','

te_matrix[cols] = te_matrix[cols].astype(float)

# filter records with MITE insertions (homo- and heterozygous, interpreted as frequency > 0.3
te_matrix_ins = te_matrix[te_matrix[cols].ge(0.3).any(axis=1)].reset_index(drop=True)
print(te_matrix_ins.to_string(max_rows=30))

# calculate sum of MITEs for every family
te_family_no = te_matrix.groupby(['family']).size()
te_family_no_sort_index = te_family_no.index.to_series().str.split('|',expand=True)
te_family_no_sort_index[1] = te_family_no_sort_index[1].astype(int)
te_family_no_sort_index = te_family_no_sort_index.sort_values([0, 1], ascending=[True, True])
te_family_no = te_family_no.reindex(te_family_no_sort_index.index).astype(int)

te_family_no_sorted = te_family_no.sort_values(ascending=False)
print(te_family_no.to_string(), type(te_family_no), te_family_no_sorted.to_string())

# te_matrix['P1-6'] = te_matrix.apply(lambda row: row['family'] if row['P1-6'] > 0 else row['P1-6'], axis=1)
# print(te_matrix.to_string(max_rows=30))

#
te_matrix[cols] = te_matrix[cols].mask(te_matrix[cols].apply(pd.to_numeric, errors='coerce').gt(0),
                         te_matrix['family'], axis=0)

# slower than the previous method
# for col in cols:
#     te_matrix[col] = te_matrix.apply(lambda row: row['family'] if row[col] > 0 else row[col], axis=1)

print(te_matrix.to_string(max_rows=30))

te_in_samples = pd.DataFrame()
for col in cols:
    te_in_sample = te_matrix.groupby([col]).size()
    te_in_samples = te_in_samples.append(te_in_sample, ignore_index=True)

te_in_samples = (te_in_samples.transpose().rename(columns={0: 'P1-6', 1: 'P1-12', 2: 'P1-22', 3: 'P1-25',
                                                          4: 'P1-26', 5: 'P1-28', 6: 'P1-88', 7: 'P1-89', 8: 'P1-90',
                                                          9: 'P1-92', 10: 'P1-93', 11: 'P1-95'}, index={0.0: 'COUNT'}).
                 fillna(0).astype(int))

print(te_in_samples.to_string())

te_in_samples_sort_index = te_in_samples.index.to_series().str.split('|',expand=True).fillna(0)
print(te_in_samples_sort_index)
te_in_samples_sort_index[1] = te_in_samples_sort_index[1].astype(int)
te_in_samples_sort_index = te_in_samples_sort_index.sort_values([0, 1], ascending=[True, True])
te_in_samples = te_in_samples.reindex(te_in_samples_sort_index.index).rename(index={'COUNT': 'SUM'})

te_in_samples = (pd.concat([te_family_no, te_in_samples], axis=1).fillna(0).astype(int).
                 rename(columns={0: 'no_of_TE_ins'}))

te_in_samples.loc['SUM', 'no_of_TE_ins'] = te_in_samples['no_of_TE_ins'].sum()
te_in_samples.columns.name = 'family'

print(te_in_samples.to_string())

te_in_samples.to_csv('P1_MITEs_count.csv',  sep='\t')

print("--- %s seconds ---" % (time.time() - start_time))
