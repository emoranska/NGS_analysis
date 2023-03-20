import pandas as pd
from ref_f_r_pairs import f_r_pairs

matrix_ref = pd.read_csv("../P1_P4/P1_F_R_pairs_merged1000.csv", sep="\t")

# creating matrix for reference insertions
for ref in f_r_pairs.itertuples():
    for x in matrix_ref.itertuples():
        if ref.chr == x.chr and ref.family == x.family and ref.pos_f in range(int(x.start), int(x.end+1)) and \
                ref.pos_r in range(int(x.start), int(x.end+1)):
            matrix_ref.at[x.Index, ref.sample] = ref.freq_average

matrix_ref.columns = matrix_ref.columns.astype("str")
matrix_ref_final = matrix_ref.reindex(columns=['chr', 'family', 'start', 'end', '1', '2', '3', '4', '5', '6', '7',
                                               '8', '9', '10', '11', '12'])
matrix_ref_final_2 = matrix_ref_final.fillna(0)
matrix_ref_final_2.to_csv('P1_ref_matrix_all.csv', sep='\t', index=False)
print(matrix_ref_final_2.to_string(max_rows=20))
