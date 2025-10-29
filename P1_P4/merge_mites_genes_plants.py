import pandas as pd

mites_matrix = pd.read_csv('../files/P4_all_MITEs_no_dupl_matrix.csv', sep='\t')
mites_in_el10 = (pd.read_csv('../files/P4_genes_with_MITEs.csv', sep='\t').
                 drop(columns=['start_cds', 'end_cds', 'start_utr', 'end_utr']))

print(mites_matrix.to_string(max_rows=30))
print(mites_in_el10.to_string(max_rows=30))

mites_in_el10_plants = (mites_in_el10.merge(mites_matrix, how="left").reset_index(drop=True))
print(mites_in_el10_plants.to_string(max_rows=100))
mites_in_el10_plants.to_csv('../files/P4_genes_with_MITEs_plants.csv', sep='\t', index=False)
