import numpy as np
import pandas as pd
import time

start_time = time.time()

# for P1
p1_mites_matrix_all = pd.read_csv('../files/P1_matrix_all.csv', sep='\t')
p1_mites_matrix_ref = pd.read_csv('../files/P1_ref_matrix_sort_2000_all.csv', sep='\t')
p1_mites_matrix_nonref = pd.read_csv('../files/P1_nonref_matrix_all.csv', sep='\t')

# for P4
p4_mites_matrix_all = pd.read_csv('../files/P4_matrix_all.csv', sep='\t')
p4_mites_matrix_ref = pd.read_csv('../files/P4_ref_matrix_sort_2000_all.csv', sep='\t')
p4_mites_matrix_nonref = pd.read_csv('../files/P4_nonref_matrix_all.csv', sep='\t')


def columns_filter(mites_matrix, pop_symbol):
    pop_cols = mites_matrix.filter(like=pop_symbol).columns
    mites_matrix[pop_cols] = mites_matrix[pop_cols].astype(float)
    return pop_cols


cols = columns_filter(p1_mites_matrix_all, 'P1-')


def mites_counts_all(te_matrix):
    # calculate sum of MITEs for every family
    te_family_no = te_matrix.groupby(['family']).size()
    te_family_no_sorted = te_family_no.sort_values(ascending=False)
    print(te_family_no.to_string(), te_family_no_sorted.to_string())

    # put family name instead of value if > 0
    te_matrix[cols] = te_matrix[cols].mask(te_matrix[cols].apply(pd.to_numeric, errors='coerce').gt(0),
                                           te_matrix['family'], axis=0)

    # slower than the previous method
    # for col in cols:
    #     te_matrix[col] = te_matrix.apply(lambda row: row['family'] if row[col] > 0 else row[col], axis=1)

    print(te_matrix.to_string(max_rows=30))

    # calculate MITEs in every sample for all families
    te_in_samples = pd.DataFrame()

    for col in cols:
        te_in_sample = te_matrix.groupby([col]).size()
        te_in_samples = pd.concat([te_in_samples, te_in_sample], axis=1)

    # for P1
    te_in_samples.columns = ['P1-6', 'P1-12', 'P1-22', 'P1-25', 'P1-26', 'P1-28', 'P1-88', 'P1-89', 'P1-90', 'P1-92',
                             'P1-93', 'P1-95']

    # for P4
    # te_in_samples.columns = ['P4-1', 'P4-7', 'P4-15', 'P4-29', 'P4-32', 'P4-35', 'P4-46', 'P4-47', 'P4-56', 'P4-60',
    #                          'P4-62', 'P4-64']

    te_in_samples = (te_in_samples.rename(index={0.0: 'COUNT'}).fillna(0).astype(int))

    print(te_in_samples.to_string())

    # sorting results
    te_in_samples_sort_index = te_in_samples.index.to_series().str.split('|', expand=True).fillna(0)
    print(te_in_samples_sort_index)
    te_in_samples_sort_index[1] = te_in_samples_sort_index[1].astype(int)

    te_in_samples_sort_index = te_in_samples_sort_index.sort_values(
        by=[0, 1], ascending=[True, True], key=lambda x: x if np.issubdtype(x.dtype, np.number) else x.str.lower())
    print(te_in_samples_sort_index)

    te_in_samples = te_in_samples.reindex(te_in_samples_sort_index.index)
    print(te_in_samples)

    te_family_no = te_family_no.reindex(te_in_samples_sort_index.index)
    print(te_family_no)

    # add column with sum of MITES
    te_in_samples = (pd.concat([te_family_no, te_in_samples], axis=1).fillna(0).astype(int).
                     rename(columns={0: 'no_of_TE_ins'}, index={'COUNT': 'SUM'}))

    te_in_samples.loc['SUM', 'no_of_TE_ins'] = te_in_samples['no_of_TE_ins'].sum()
    te_in_samples.columns.name = 'family'

    print(te_in_samples.to_string())

    # saving as .csv
    # te_in_samples.to_csv('../files/P4_MITEs_nonref_count.csv', sep='\t')


# mites_counts_all(p1_mites_matrix_all)


def mites_counts_ins(te_matrix):

    # for homozygous insertions - put family name instead of value if >= 0.7
    # te_matrix[cols] = te_matrix[cols].mask(te_matrix[cols].apply(pd.to_numeric, errors='coerce').ge(0.7),
    #                                        te_matrix['family'], axis=0)

    # for heterozygous insertions - put family name instead of value if  < 0.3 >= 0.7
    for col in cols:
        te_matrix[col] = te_matrix.apply(lambda row: row['family'] if 0.3 <= row[col] < 0.7 else row[col], axis=1)

    print(te_matrix.to_string(max_rows=30))

    # calculate MITEs in every sample for all families
    te_in_samples = pd.DataFrame()

    for col in cols:
        te_in_sample = te_matrix.groupby([col]).size()
        te_in_samples = pd.concat([te_in_samples, te_in_sample], axis=1)

    # for P1
    te_in_samples.columns = ['P1-6', 'P1-12', 'P1-22', 'P1-25', 'P1-26', 'P1-28', 'P1-88', 'P1-89', 'P1-90', 'P1-92',
                             'P1-93', 'P1-95']

    te_in_samples = te_in_samples.fillna(0)
    te_in_samples = te_in_samples.filter(regex='^\D', axis=0).astype(int)

    print(te_in_samples.to_string())

    # sorting results
    te_in_samples_sort_index = te_in_samples.index.to_series().str.split('|', expand=True).fillna(0)
    print(te_in_samples_sort_index)
    te_in_samples_sort_index[1] = te_in_samples_sort_index[1].astype(int)

    te_in_samples_sort_index = te_in_samples_sort_index.sort_values(
        by=[0, 1], ascending=[True, True], key=lambda x: x if np.issubdtype(x.dtype, np.number) else x.str.lower())
    print(te_in_samples_sort_index)

    te_in_samples = te_in_samples.reindex(te_in_samples_sort_index.index)
    print(te_in_samples.to_string())

    # saving as .csv
    te_in_samples.to_csv('../files/P1_MITEs_heterozyg_ins_count.csv', sep='\t')


mites_counts_ins(p1_mites_matrix_all)

print("--- %s seconds ---" % (time.time() - start_time))
