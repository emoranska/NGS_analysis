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

te_matrix[cols] = te_matrix[cols].astype(float)

# create column with list of samples with MITE insertions, interpreted as frequency > 0.7
#zrobić odwrotność i będzie super
te_matrix_ins = te_matrix[te_matrix[cols].ge(0.9).any(axis=1)]
print(te_matrix_ins.to_string(max_rows=30))
