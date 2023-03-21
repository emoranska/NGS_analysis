import pandas as pd
import time

start_time = time.time()

matrix_nonref = pd.read_csv('P1_nonref_merged1000.csv', sep='\t')
nonref_all_sort = pd.read_csv('P1_nonref_all.csv', sep='\t')

# creating matrix for nonref TE insertions
for te in nonref_all_sort.itertuples():
    for x in matrix_nonref.itertuples():
        if te.chr == x.chr and te.family == x.family and te.pos in range(int(x.start), int(x.end+1)):
            matrix_nonref.at[x.Index, te.sample] = te.frequency

matrix_nonref.columns = matrix_nonref.columns.astype("str")
matrix_nonref_final = matrix_nonref.reindex(columns=['chr', 'family', 'start', 'end', '1', '2', '3', '4', '5', '6',
                                                     '7', '8', '9', '10', '11', '12']).fillna(0)

print(matrix_nonref_final.to_string(max_rows=20))
matrix_nonref_final.to_csv('P1_nonref_matrix_all.csv', sep='\t', index=False)

print("--- %s seconds ---" % (time.time() - start_time))
