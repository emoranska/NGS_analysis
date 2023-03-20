import pandas as pd

# concat all parts of unpaired nonref F and R results with FR file
fr = pd.read_csv('../P1_P4/P4_FR.csv', sep='\t')
un_chr1 = pd.read_csv('P4_nonref_unpaired_chr/P4_nonref_unpaired_chr1.csv', sep="\t")
un_chr2 = pd.read_csv('P4_nonref_unpaired_chr/P4_nonref_unpaired_chr2.csv', sep="\t")
un_chr3 = pd.read_csv('P4_nonref_unpaired_chr/P4_nonref_unpaired_chr3.csv', sep="\t")
un_chr4 = pd.read_csv('P4_nonref_unpaired_chr/P4_nonref_unpaired_chr4.csv', sep="\t")
un_chr5 = pd.read_csv('P4_nonref_unpaired_chr/P4_nonref_unpaired_chr5.csv', sep="\t")
un_chr6 = pd.read_csv('P4_nonref_unpaired_chr/P4_nonref_unpaired_chr6.csv', sep="\t")
un_chr7 = pd.read_csv('P4_nonref_unpaired_chr/P4_nonref_unpaired_chr7.csv', sep="\t")
un_chr8 = pd.read_csv('P4_nonref_unpaired_chr/P4_nonref_unpaired_chr8.csv', sep="\t")
un_chr9 = pd.read_csv('P4_nonref_unpaired_chr/P4_nonref_unpaired_chr9.csv', sep="\t")

nonref_final = pd.concat([fr, un_chr1, un_chr2, un_chr3, un_chr4, un_chr5, un_chr6, un_chr7, un_chr8, un_chr9],
                         ignore_index=True).sort_values(by=['chr', 'pos'])

nonref_final.to_csv('P4_nonref_all.csv', index=False, sep='\t')
print(nonref_final.to_string(max_rows=20))
