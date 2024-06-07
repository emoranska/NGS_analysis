import pandas as pd
import time
from pathlib import Path
import ast

start_time = time.time()

te_matrix = pd.read_csv('../files/P1_matrix_all.csv', sep='\t')
# print(te_matrix.to_string(max_rows=30))

cols = te_matrix.filter(like='P1-').columns
sep = ','

# create column with list of samples with MITE insertions, interpreted as frequency > 0.8
te_matrix['ok_ins'] = (pd.Series(te_matrix[cols].ge(0.8).dot((cols + sep))).str.rstrip(sep).str.split(',').
                       sort_values().apply(lambda x: sorted(x)))

# create column with list of samples without MITE insertions, interpreted as frequency < 0.3
te_matrix['less_0.3'] = (pd.Series(te_matrix[cols].le(0.3).dot((cols + sep))).str.rstrip(sep).str.split(',').
                         sort_values().apply(lambda x: sorted(x)))
'''
# create column with list of samples with frequency == 0
te_matrix['eq_0'] = (pd.Series(te_matrix[cols].eq(0).dot((cols + sep))).str.rstrip(sep).str.split(',').sort_values().
                     apply(lambda x: sorted(x)))

# create column with list of samples without MITE insertions --> frequency < 0.3 and != 0
notna_only = te_matrix[['family']].notna().all(axis=1)
te_matrix.loc[notna_only, 'no_ins'] = [list(sorted(set(a) - set(b))) for a, b in
                                       zip(te_matrix.loc[notna_only, 'less_0.3'], te_matrix.loc[notna_only, 'eq_0'])]
'''

# create columns with list of 3 first samples with and without insertions
te_matrix['ins_3first'] = te_matrix['ok_ins'].apply(lambda y: y[0:3])
te_matrix['no_ins_3first'] = te_matrix['less_0.3'].apply(lambda i: i[0:3])

print(te_matrix.to_string(max_rows=30))

te_filter = te_matrix.drop(columns=['less_0.3', 'ok_ins']).rename(columns={
    'ins_3first': 'te_ins', 'no_ins_3first': 'no_te_ins'})

print(te_filter.to_string(max_rows=50))

'''
# checking data type in new columns because of problems with empty strings inside lists
print(te_filter.loc[16, 'no_te_ins'], type(te_filter.loc[16, 'no_te_ins']), len(te_filter.loc[16, 'no_te_ins']))
for x in te_filter.loc[16, 'no_te_ins']:
    print(x, type(x), len(x))

print(te_filter.loc[13, 'no_te_ins'], type(te_filter.loc[13, 'no_te_ins']), len(te_filter.loc[13, 'no_te_ins']))
for x in te_filter.loc[13, 'no_te_ins']:
    print(x, type(x), len(x), x.count(''))

print(te_filter.loc[1, 'te_ins'], type(te_filter.loc[1, 'te_ins']), len(te_filter.loc[1, 'te_ins']))
for x in te_filter.loc[1, 'te_ins']:
    print(x, type(x), len(x), x.count(''))

print(te_filter.loc[13, 'te_ins'], type(te_filter.loc[13, 'te_ins']), len(te_filter.loc[13, 'te_ins']))
for x in te_filter.loc[13, 'te_ins']:
    print(x, len(x), x.count(''))
    print(type(x))
'''

# filter out empty strings from list of samples in 'te_ins' columns !!!
te_filter['te_ins'] = te_filter['te_ins'].apply(lambda x: list(filter(None, x)))
te_filter['no_te_ins'] = te_filter['no_te_ins'].apply(lambda x: list(filter(None, x)))

# checking previous line effect
for x in te_filter.loc[13, 'te_ins']:
    print(x, len(x), x.count(''))
    print(type(x))

# first method for filtering only rows containing 3 samples with MITE insertions AND 3 samples without it
te_filter_final = te_filter[(te_filter['te_ins'].map(lambda d: len(d) == 3)) &
                            (te_filter['no_te_ins'].map(lambda d: len(d) == 3))].reset_index(drop=True)

# second method, slightly slower
# te_filter_final = te_filter[(te_filter['te_ins'].str.len() == 3) & (te_filter['no_te_ins'].str.len() == 3)]

te_filter_final = te_filter_final.drop(columns=te_filter_final.filter(like='P1-').columns).reset_index(drop=True)
print(te_filter_final.to_string(max_rows=30))

# te_filter_final.to_csv('../files/P1_MITE_list.csv', sep='\t', index=False)

print("--- %s seconds ---" % (time.time() - start_time))
