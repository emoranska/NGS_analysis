import pandas as pd
import numpy as np

bins_and_de = (pd.read_csv('../files/P1_all_EL10_genes_in_all_bins_with_DE_and_kallisto.csv', sep='\t').
               drop(columns=['start_gene', 'end_gene'])).rename(columns={'Start': 'start_gene', 'End': 'end_gene'})
# print(bins_and_de.memory_usage(deep=True))
start_gene = bins_and_de.pop('start_gene')
bins_and_de.insert(5, start_gene.name, start_gene)

end_gene = bins_and_de.pop('end_gene')
bins_and_de.insert(6, end_gene.name, end_gene)

print(bins_and_de.to_string(max_rows=30))

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

mites_with_unique_no = (((opt_mites_all.merge(sets_list, how='cross').
                        query('((te_ins == one_one) & (no_te_ins == zero_zero)) | '
                              '((te_ins == zero_zero) & (no_te_ins == one_one))')).
                        drop(columns=['one_one', 'zero_zero']).rename(columns={'start': 'start_te', 'end': 'end_te'})).
                        reset_index(drop=True))

print(mites_with_unique_no.to_string(max_rows=30))


mites_in_genes = (bins_to_compare.merge(mites_with_unique_no, how='outer', on=['chr', 'unique_no']).
                  query('(start_te >= start_gene) & (end_te <= end_gene)').reset_index(drop=True))
mites_in_genes['mite_loc'] = 'intron'
print('MITEs in genes:', '\n', mites_in_genes.to_string(max_rows=30))

exons_el10 = pd.read_csv('../files/EL10_exons_bed.csv', sep='\t').drop(columns=['gene_name'])
mites_in_exons = (mites_in_genes.merge(exons_el10, how='outer', on=['chr']).
                  query('(start_te >= start_ex & end_te <= end_ex)').drop_duplicates().reset_index(drop=True))

print('MITEs in exons:', '\n', mites_in_exons.to_string(max_rows=30))

mites_in_genes_ex = mites_in_genes.merge(mites_in_exons, how='outer')
mites_in_genes_ex.loc[~np.isnan(mites_in_genes_ex['start_ex']), 'mite_loc'] = 'exon'
# df.loc[df['salary'] > 50,'is_rich_method3'] = 'yes'
print('MITEs in genes with exons checked:', '\n', mites_in_genes_ex.to_string())

mites_on_board = (bins_to_compare.merge(mites_with_unique_no, how='outer', on=['chr', 'unique_no']).
                  query('(start_te < start_gene) & (end_te >= start_gene) | '
                        '(end_te > end_gene) & (start_te <= end_gene)').
                  reset_index(drop=True))
print('MITEs on board:', '\n', mites_on_board.to_string(max_rows=30))

mites_ob_in_exons = (mites_on_board.merge(exons_el10, how='outer', on=['chr']).
                     query('(start_te < start_ex & end_te > start_ex) | '
                            '(start_te < end_ex & end_te > end_ex)').drop_duplicates().reset_index(drop=True))
print('MITEs on board in exons:', '\n', mites_ob_in_exons.to_string(max_rows=30))
# df['is_rich_method2'] = ['yes' if x >= 50 else 'no' for x in df['salary']]
mites_on_board_exons = mites_on_board.merge(mites_ob_in_exons, how='outer')

for x, y, z in zip(mites_on_board_exons['start_ex'], mites_on_board_exons['start_te'], mites_on_board_exons['gene_strand']):
    if np.isnan(x) and ((y < x and z == '+') or (y > x and z == '-')):
        mites_on_board_exons['mite_loc'] = 'board_up'
    elif np.isnan(x) and ((y < x and z == '-') or (y > x and z == '+')):
        mites_on_board_exons['mite_loc'] = 'board_down'
    elif ~np.isnan(x) and ((y < x and z == '+') or (y > x and z == '-')):
        mites_on_board_exons['mite_loc'] = 'exon_up'
    elif ~np.isnan(x) and ((y < x and z == '-') or (y > x and z == '+')):
        mites_on_board_exons['mite_loc'] = 'exon_down'

print(mites_on_board_exons.to_string())

mites_updown = (bins_to_compare.merge(mites_with_unique_no, how='outer', on=['chr', 'unique_no']).
                query('(start_te >= start_gene - 2000) & (end_te < start_gene) | (end_te <= end_gene + 2000) & '
                      '(start_te > end_gene)').
                reset_index(drop=True))
print('MITEs updown:', '\n', mites_updown.to_string(max_rows=30))


bins_genes_mites = (opt_bins_and_de.merge(mites_in_genes, how='outer', on=['chr', 'start_gene', 'end_gene', 'unique_no',
                                                                  'set_no']).reset_index(drop=True))
# print(bins_genes_mites.to_string(max_rows=200))

te_ref = pd.read_csv('../files/P1_ref_matrix_sort_2000_all.csv', sep='\t')
te_ref = te_ref.drop(columns=te_ref.columns[-12:], axis=1)
# print(te_ref.to_string(max_rows=30))
te_nonref = pd.read_csv('../files/P1_nonref_matrix_all.csv', sep='\t')
te_nonref = te_nonref.drop(columns=te_nonref.columns[-12:], axis=1)
# print(te_nonref.to_string(max_rows=30))

ref_check = bins_genes_mites.merge(te_ref, on=['chr', 'family', 'start_te', 'end_te'])
ref_check['te_type'] = 'ref'
print('Ref check:', '\n', ref_check.to_string(max_rows=30))

nonref_check = bins_genes_mites.merge(te_nonref, on=['chr', 'family', 'start_te', 'end_te'])
nonref_check['te_type'] = 'nonref'
print('Nonref check:', '\n', nonref_check.to_string(max_rows=30))

type_checked = pd.concat([ref_check, nonref_check]).reset_index(drop=True)
print('Type checked:', '\n', type_checked.to_string(max_rows=30))

# mites_in_genes_te_types = mites_in_genes.merge(type_checked, how='outer')
# print(mites_in_genes_te_types.to_string(max_rows=100))


type_checked_exons = (type_checked.merge(exons_el10, how='outer', on=['chr']).
                      query('(start_te >= start_ex & end_te <= end_ex) | (start_te < start_ex & end_te > start_ex) | '
                            '(start_te < end_ex & end_te > end_ex)')).reset_index(drop=True)
print(type_checked_exons.to_string(max_rows=30))

type_checked_ex_int = type_checked.merge(type_checked_exons, how='outer').drop_duplicates().reset_index(drop=True)
print(type_checked_ex_int.to_string(max_rows=30))

# print(mites_in_genes.to_string(max_rows=50))
# all_together.to_csv('../files/P1_alltogether.csv', sep='\t', index=False)
