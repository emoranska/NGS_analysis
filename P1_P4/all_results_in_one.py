import pandas as pd

bins_and_de = pd.read_csv('../files/P1_all_EL10_genes_in_all_bins_with_DE_and_kallisto.csv', sep='\t')
# print(bins_and_de.memory_usage(deep=True))

bins_and_de_int = bins_and_de.select_dtypes(include=['int'])
bins_conv_int = bins_and_de_int.apply(pd.to_numeric, downcast='unsigned')
# print(bins_and_de_int.memory_usage(deep=True))
# print(bins_conv_int.memory_usage(deep=True))

bins_and_de_float = bins_and_de.select_dtypes(include=['float'])
bins_conv_float = bins_and_de_float.apply(pd.to_numeric, downcast='float')
# print(bins_and_de_float.memory_usage(deep=True))
# print(bins_conv_float.memory_usage(deep=True))

bins_and_de_obj = bins_and_de.select_dtypes(include=['object'])
bins_conv_obj = pd.DataFrame()
for col in bins_and_de_obj.columns:
    num_unique_values = len(bins_and_de_obj[col].unique())
    num_total_values = len(bins_and_de_obj[col])
    if num_unique_values / num_total_values < 0.5:
        bins_conv_obj.loc[:, col] = bins_and_de_obj[col].astype('category')
    else:
        bins_conv_obj.loc[:, col] = bins_and_de_obj[col]

opt_bins_and_de = bins_and_de.copy()
opt_bins_and_de[bins_conv_int.columns] = bins_conv_int
opt_bins_and_de[bins_conv_float.columns] = bins_conv_float
opt_bins_and_de[bins_conv_obj.columns] = bins_conv_obj
print(bins_and_de.memory_usage(deep=True).sum())
print('Memory usage: {:03.2f}'.format(bins_and_de.memory_usage(deep=True).sum()))
print(opt_bins_and_de.memory_usage(deep=True).sum())
print(bins_conv_obj.memory_usage(deep=True).sum())
print(opt_bins_and_de.to_string(max_rows=30))
print(bins_conv_obj.to_string(max_rows=30))

bins_to_compare = opt_bins_and_de[['chr', 'start_gene', 'end_gene', 'unique_no']]
# bins_to_compare['bins_no'] = bins_to_compare.index
print(bins_to_compare.to_string(max_rows=30))
print(bins_to_compare.memory_usage(deep=True).sum())


# mites_padj05 = pd.read_csv('../files/P1_correct_MITES_in_genes_and_updown.csv', sep='\t',
#                            dtype={'log2FoldChange': 'float64', 'lfcSE': 'float64'})
# print(mites_padj05.to_string(max_rows=30))


samples = pd.read_csv('../files/P1_list_to_DE_20more_genes.csv', sep='\t')
sets_list = (pd.read_csv('../files/P1_sets_test.csv', sep='\t').
             sort_values(by=['set_no', 'unique_no']).drop(columns=['genes_count']))
print(sets_list.to_string())

mites_all = pd.read_csv('../files/P1_MITEs_with_bins_sets.csv', sep='\t')

mites_all_int = mites_all.select_dtypes(include=['int'])
mites_conv_int = mites_all_int.apply(pd.to_numeric, downcast='unsigned')
print(mites_all_int.memory_usage(deep=True).sum())
print(mites_conv_int.memory_usage(deep=True).sum())

mites_all_float = mites_all.select_dtypes(include=['float'])
mites_conv_float = mites_all_float.apply(pd.to_numeric, downcast='float')
print(mites_all_float.memory_usage(deep=True))
print(mites_conv_float.memory_usage(deep=True))

opt_mites_all = mites_all.copy()
opt_mites_all[mites_conv_int.columns] = mites_all_int
opt_mites_all[mites_conv_float.columns] = mites_all_float

mites_all_obj = mites_all.select_dtypes(include=['object']).copy()
mites_conv_obj = pd.DataFrame()
for col in mites_all_obj.columns:
    num_unique_values = len(mites_all_obj[col].unique())
    num_total_values = len(mites_all_obj[col])
    if num_unique_values / num_total_values < 0.5:
        mites_conv_obj.loc[:, col] = mites_all_obj[col].astype('category')
    else:
        mites_conv_obj.loc[:, col] = mites_all_obj[col]

opt_mites_all[mites_conv_obj.columns] = mites_conv_obj
print((mites_all.memory_usage(deep=True).sum()))
print((opt_mites_all.memory_usage(deep=True).sum()))
print(mites_conv_obj.memory_usage(deep=True).sum())
print(opt_mites_all.to_string(max_rows=30))

# mites_conv_obj['te_no'] = mites_conv_obj.index
# mites_conv_obj['te_no'] = mites_conv_obj['te_no'].apply(pd.to_numeric, downcast='unsigned')
print(mites_conv_obj.memory_usage(deep=True).sum())
print(mites_conv_obj.to_string(max_rows=30))

mites_with_unique_no = (opt_mites_all.merge(sets_list, how='cross').
                        query('((te_ins == one_one) & (no_te_ins == zero_zero)) | '
                              '((te_ins == zero_zero) & (no_te_ins == one_one))')).drop(columns=['one_one', 'zero_zero'])

print(mites_with_unique_no.to_string(max_rows=30))

part_df = (bins_to_compare.merge(mites_with_unique_no, how='outer', on=['chr', 'unique_no']).
           query('(start >= start_gene) & (end <= end_gene)').reset_index(drop=True))
print(part_df.to_string(max_rows=30))

mites_in_genes = (opt_bins_and_de.merge(part_df, how='outer', on=['chr', 'start_gene', 'end_gene', 'unique_no', 'set_no']).
                  rename(columns={'start': 'start_te', 'end': 'end_te'}).reset_index(drop=True))
print(mites_in_genes.to_string(max_rows=300))

te_ref = pd.read_csv('../files/P1_ref_matrix_sort_2000_all.csv', sep='\t')
te_ref = te_ref.drop(columns=te_ref.columns[-12:], axis=1)
print(te_ref.to_string(max_rows=30))
te_nonref = pd.read_csv('../files/P1_nonref_matrix_all.csv', sep='\t')
te_nonref = te_nonref.drop(columns=te_nonref.columns[-12:], axis=1)
print(te_nonref.to_string(max_rows=30))

ref_check = mites_in_genes.merge(te_ref, on=['chr', 'family', 'start_te', 'end_te'])
ref_check['te_type'] = 'ref'
print(ref_check.to_string(max_rows=30))

nonref_check = mites_in_genes.merge(te_nonref, on=['chr', 'family', 'start_te', 'end_te'])
nonref_check['te_type'] = 'nonref'
print(nonref_check.to_string(max_rows=30))
# mites_in_genes =

# print(mites_in_genes.to_string(max_rows=50))
# all_together.to_csv('../files/P1_alltogether.csv', sep='\t', index=False)
