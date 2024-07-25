import ast

import pandas as pd
import time
import ast

start_time = time.time()

sample_sets = pd.read_csv('../files/P1_list_to_DE_20more_genes.csv', sep='\t')
sample_sets['one_one'] = sample_sets['one_one'].apply(lambda x: ast.literal_eval(x))
sample_sets['zero_zero'] = sample_sets['zero_zero'].apply(lambda x: ast.literal_eval(x))
print(sample_sets.iloc[0,2], type(sample_sets.iloc[0,2]))
print(sample_sets.to_string(max_rows=20))

te_matrix = pd.read_csv('../files/P1_matrix_all.csv', sep='\t')
# te_matrix = pd.read_csv('../files/P4_matrix_all.csv', sep='\t')
# print(te_matrix.to_string(max_rows=30))

te_matrix['start'] = te_matrix['start'].astype(int)
te_matrix['end'] = te_matrix['end'].astype(int)

cols = te_matrix.filter(like='P1-').columns
# cols = te_matrix.filter(like='P4-').columns
sep = ','

# create column with list of samples with homozygous MITE insertions, interpreted as frequency > 0.7
te_matrix['ok_ins'] = (pd.Series(te_matrix[cols].ge(0.7).dot((cols + sep))).str.rstrip(sep).str.split(',').
                       sort_values().apply(lambda t: sorted(t)))

# create column with list of samples without MITE insertions (homozygous), interpreted as frequency < 0.3
te_matrix['no_ins'] = (pd.Series(te_matrix[cols].le(0.3).dot((cols + sep))).str.rstrip(sep).str.split(',').
                         sort_values().apply(lambda z: sorted(z)))

# filter out empty strings from list of samples in 'ok_ins' and 'no_ins' columns !!!
te_matrix['ok_ins'] = te_matrix['ok_ins'].apply(lambda x: list(filter(None, x)))
te_matrix['no_ins'] = te_matrix['no_ins'].apply(lambda x: list(filter(None, x)))

te_matrix = te_matrix[(te_matrix['ok_ins'].map(lambda d: len(d) >= 3)) &
                             (te_matrix['no_ins'].map(lambda d: len(d) >= 3))].reset_index(drop=True)

print(te_matrix.iloc[0, 17], type(te_matrix.iloc[0, 17]))
print(te_matrix.to_string(max_rows=30))

matches = [x for x in sample_sets.iloc[0,2] if x in te_matrix.iloc[0, 17]]
print(f'Matches: {matches}')

sample_sets = sample_sets.drop(columns=['df_unique_index', 'genes_count'])
out = te_matrix.merge(sample_sets, how='cross')
print(out.to_string(max_rows=50))

cols_s = list(sample_sets)
print(cols_s)

ins = out['ok_ins'].apply(set)
print(ins)

for c in cols_s:
    out[f'{c}_ins'] = [[x for x in lst if x in ref] for ref, lst in zip(ins, out[c])]

print(out.to_string(max_rows=50))

no_ins = out['no_ins'].apply(set)

for c in cols_s:
    out[f'{c}_no_ins'] = [[x for x in lst if x in ref] for ref, lst in zip(no_ins, out[c])]

print(out.to_string(max_rows=50))

mites_in_bins = out[((out['one_one_ins'].map(lambda d: len(d) == 3)) & (out['zero_zero_no_ins'].map(lambda d:
                                                                                                    len(d) == 3))) |
                    ((out['zero_zero_ins'].map(lambda d: len(d) == 3)) & (out['one_one_no_ins'].map(lambda d:
                                                                                                    len(d) == 3)))]

mites_in_bins = (mites_in_bins.drop(columns=cols).
                 drop(columns=['ok_ins', 'no_ins']).reset_index(drop=True))

mites_in_bins['ins'] = mites_in_bins['one_one_ins'] + mites_in_bins['zero_zero_ins']
mites_in_bins['no_ins'] = mites_in_bins['one_one_no_ins'] + mites_in_bins['zero_zero_no_ins']

print(mites_in_bins.to_string(max_rows=50))

mites_in_bins_final = (mites_in_bins.drop(columns=['one_one', 'zero_zero', 'one_one_ins', 'zero_zero_ins',
                                                  'one_one_no_ins', 'zero_zero_no_ins']).rename
                       (columns={'ins': 'te_ins', 'no_ins': 'no_te_ins'}))
print(mites_in_bins_final.to_string(max_rows=50))
mites_in_bins_final.to_csv('../files/P1_MITEs_in_bins.csv', sep='\t')

'''
# create columns with list of 3 first samples with and without insertions - WRONG! Simplified too much!
te_matrix['ins_3first'] = te_matrix['ok_ins'].apply(lambda y: y[0:3])
te_matrix['no_ins_3first'] = te_matrix['no_ins'].apply(lambda i: i[0:3])

# print(te_matrix.to_string(max_rows=30))

te_filter = te_matrix.drop(columns=['no_ins', 'ok_ins']).rename(columns={
    'ins_3first': 'te_ins', 'no_ins_3first': 'no_te_ins'})

# print(te_filter.to_string(max_rows=50))

# filter out empty strings from list of samples in 'te_ins' and 'no_te_ins' columns !!!
te_filter['te_ins'] = te_filter['te_ins'].apply(lambda x: list(filter(None, x)))
te_filter['no_te_ins'] = te_filter['no_te_ins'].apply(lambda x: list(filter(None, x)))

# first method for filtering only rows containing 3 samples with MITE insertions AND 3 samples without it
te_filter_final = te_filter[(te_filter['te_ins'].map(lambda d: len(d) == 3)) &
                            (te_filter['no_te_ins'].map(lambda d: len(d) == 3))].reset_index(drop=True)

# second method, slightly slower
# te_filter_final = te_filter[(te_filter['te_ins'].str.len() == 3) & (te_filter['no_te_ins'].str.len() == 3)]

te_filter_final = te_filter_final.drop(columns=te_filter_final.filter(like='P1-').columns).reset_index(drop=True)
# te_filter_final = te_filter_final.drop(columns=te_filter_final.filter(like='P4-').columns).reset_index(drop=True)
# print(te_filter_final.to_string(max_rows=30))

# te_filter_final.to_csv('../files/P1_MITE_list.csv', sep='\t', index=False)
# te_filter_final.to_csv('../files/P4_MITE_list.csv', sep='\t', index=False)
'''

print("--- %s seconds ---" % (time.time() - start_time))
