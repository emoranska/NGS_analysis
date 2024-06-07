import pandas as pd
import time

start_time = time.time()

te_matrix = pd.read_csv('../files/P1_matrix_all.csv', sep='\t')
# te_matrix = pd.read_csv('../files/P4_matrix_all.csv', sep='\t')
# print(te_matrix.to_string(max_rows=30))

te_matrix['start'] = te_matrix['start'].astype(int)
te_matrix['end'] = te_matrix['end'].astype(int)

cols = te_matrix.filter(like='P1-').columns
# cols = te_matrix.filter(like='P4-').columns
sep = ','

# create column with list of samples with MITE insertions, interpreted as frequency > 0.8
te_matrix['ok_ins'] = (pd.Series(te_matrix[cols].ge(0.7).dot((cols + sep))).str.rstrip(sep).str.split(',').
                       sort_values().apply(lambda t: sorted(t)))

# create column with list of samples without MITE insertions, interpreted as frequency < 0.3
te_matrix['less_0.3'] = (pd.Series(te_matrix[cols].le(0.3).dot((cols + sep))).str.rstrip(sep).str.split(',').
                         sort_values().apply(lambda z: sorted(z)))

# create columns with list of 3 first samples with and without insertions
te_matrix['ins_3first'] = te_matrix['ok_ins'].apply(lambda y: y[0:3])
te_matrix['no_ins_3first'] = te_matrix['less_0.3'].apply(lambda i: i[0:3])

print(te_matrix.to_string(max_rows=30))

te_filter = te_matrix.drop(columns=['less_0.3', 'ok_ins']).rename(columns={
    'ins_3first': 'te_ins', 'no_ins_3first': 'no_te_ins'})

print(te_filter.to_string(max_rows=50))

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
print(te_filter_final.to_string(max_rows=30))

te_filter_final.to_csv('../files/P1_MITE_list.csv', sep='\t', index=False)
# te_filter_final.to_csv('../files/P4_MITE_list.csv', sep='\t', index=False)

print("--- %s seconds ---" % (time.time() - start_time))
