import pandas as pd
import numpy as np
from concat_bins_with_genes_and_rest import find_unique_no

bins_and_de = (pd.read_csv('../files/P4_all_EL10_genes_in_all_bins_with_DE_and_kallisto_tpm.csv', sep='\t'))
# # print(bins_and_de.memory_usage(deep=True))

# change gene localisation according to data from DE (Start, End)
a = bins_and_de['start_gene']
b = bins_and_de['GeneID']
c = bins_and_de['Start']
d = bins_and_de['end_gene']
e = bins_and_de['End']

cond_start_end = [pd.isnull(b), pd.notnull(b)]
choices_start = [a, c]
bins_and_de['start_test'] = np.select(cond_start_end, choices_start, default=0)
choices_end = [d, e]
bins_and_de['end_test'] = np.select(cond_start_end, choices_end, default=0)

bins_and_de = (bins_and_de.drop(columns=['start_gene', 'end_gene', 'Start', 'End']).
               rename(columns={'start_test': 'start_gene', 'end_test': 'end_gene'}))

start_gene = bins_and_de.pop('start_gene')
bins_and_de.insert(5, start_gene.name, start_gene)
bins_and_de['start_gene'] = bins_and_de['start_gene'].astype(int)

end_gene = bins_and_de.pop('end_gene')
bins_and_de.insert(6, end_gene.name, end_gene)
bins_and_de['end_gene'] = bins_and_de['end_gene'].astype(int)

print("Bins and DE results:", '\n', bins_and_de.to_string(max_rows=50))

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

bins_to_compare = opt_bins_and_de[['chr', 'start_gene', 'end_gene', 'gene_name', 'gene_strand', 'unique_no']]
print(bins_to_compare.to_string(max_rows=30))
print(bins_to_compare.memory_usage(deep=True).sum())


samples = pd.read_csv('../files/P4_list_to_DE_20more_genes.csv', sep='\t')
sets_list = (pd.read_csv('../files/P4_sets_test.csv', sep='\t').
             sort_values(by=['set_no', 'unique_no']).drop(columns=['genes_count']))
# sets_list = find_unique_no(samples).drop(columns=['genes_count'])
print(sets_list.to_string())

mites_all = pd.read_csv('../files/P4_MITEs_with_bins_sets.csv', sep='\t').drop(columns=['set_no'])

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
print(mites_conv_obj.to_string(max_rows=30))

mites_with_unique_no = (((opt_mites_all.merge(sets_list, how='cross').
                        query('((te_ins == one_one) & (no_te_ins == zero_zero)) | '
                              '((te_ins == zero_zero) & (no_te_ins == one_one))')).
                        drop(columns=['one_one', 'zero_zero']).rename(columns={'start': 'start_te', 'end': 'end_te'})).
                        reset_index(drop=True))

print('MITEs with unique_no:', '\n', mites_with_unique_no.to_string(max_rows=30))

# add MITEs to genes
mites_in_genes = (bins_to_compare.merge(mites_with_unique_no, how='outer', on=['chr', 'unique_no']).
                  query('(start_te >= start_gene) & (end_te <= end_gene)').reset_index(drop=True))
mites_in_genes['mite_loc'] = 'intron'
print('MITEs in genes:', '\n', mites_in_genes.to_string(max_rows=30))

# check MITEs in exons
exons_el10 = pd.read_csv('../files/EL10_exons_bed.csv', sep='\t').drop(columns=['gene_name', 'gene_strand'])
mites_in_exons = (mites_in_genes.merge(exons_el10, how='outer', on=['chr']).
                  query('(start_te >= start_ex & end_te <= end_ex)').drop_duplicates().reset_index(drop=True))
mites_in_exons['mite_loc'] = 'exon'
print('MITEs in exons:', '\n', mites_in_exons.to_string(max_rows=30))

cds_el10 = pd.read_csv('../files/EL10_CDS_bed.csv', sep='\t')

mites_in_exons_with_cds_check = (mites_in_exons.merge(cds_el10, how='outer', on=['chr', 'gene_name', 'gene_strand']).
                                 dropna(subset=['unique_no']).reset_index(drop=True))
print('MITEs with cds coordinants:', '\n', mites_in_exons_with_cds_check.to_string(max_rows=50))

# change 'mite_loc' for 'exon' into '5_UTR' or '3_UTR'
ste = mites_in_exons_with_cds_check['start_te']
sc = mites_in_exons_with_cds_check['start_cds']
ec = mites_in_exons_with_cds_check['end_cds']
gs = mites_in_exons_with_cds_check['gene_strand']

cond_utr_exons = [((ste.lt(sc) & gs.eq('+')) | (ste.gt(sc) & gs.eq('-'))),
                  ((ste.ge(ec) & gs.eq('+')) | (ste.lt(ec) & gs.eq('-')))]
choices_utr = ['5_UTR', '3_UTR']
mites_in_exons_with_cds_check['mite_loc'] = np.select(cond_utr_exons, choices_utr, default='exon')
mites_in_exons_with_cds_check = ((mites_in_exons_with_cds_check.
                                 drop_duplicates(subset=['chr', 'start_gene', 'end_gene', 'gene_strand', 'unique_no',
                                                         'family', 'start_te', 'end_te', 'start_ex', 'end_ex']).
                                 dropna(subset=['unique_no']).reset_index(drop=True)).
                                 rename(columns={'start_cds': 'start_annot_cds', 'end_cds': 'end_annot_cds'}))
print('MITEs in exons with UTRs checked:', '\n', mites_in_exons_with_cds_check.to_string(max_rows=200))

# dokończyć sprawdzanie UTR!!! --> potem dla board --> potem 'exon' merge with lncRNA
exon_no_cds_check = mites_in_exons_with_cds_check[mites_in_exons_with_cds_check['mite_loc'].eq('exon')]
print('MITEs in exons not in cds:', '\n', exon_no_cds_check.to_string())

mites_in_genes_ex = mites_in_genes.merge(mites_in_exons_with_cds_check, how='outer').drop_duplicates()
# mites_in_genes_ex.loc[~np.isnan(mites_in_genes_ex['start_ex']), 'mite_loc'] = 'exon'

print('MITEs in genes with exons, UTRs and cds checked:', '\n', mites_in_genes_ex.to_string(max_rows=50))

# check MITEs in cds
cds_el10_drop = cds_el10.drop(columns=['gene_name', 'gene_strand'])
mites_in_cds = (mites_in_exons.merge(cds_el10_drop, how='outer', on=['chr']).
                  query('(start_te >= start_cds & end_te <= end_cds)').drop_duplicates().reset_index(drop=True))
mites_in_cds['mite_loc'] = 'cds'
print('MITEs in cds:', '\n', mites_in_cds.to_string(max_rows=30))

mites_in_genes_ex_cds = ((mites_in_genes_ex.merge(mites_in_cds, how='outer').
                         sort_values(by=['chr', 'start_gene', 'start_te', 'mite_loc']).
                         drop_duplicates(subset=['chr', 'start_gene', 'end_gene', 'gene_strand', 'unique_no', 'family',
                                                 'start_te', 'start_ex', 'end_ex'], keep='last')).
                         drop_duplicates(subset=['chr', 'start_gene', 'end_gene', 'gene_strand', 'unique_no', 'family',
                                               'start_te', 'end_te'], keep='first').
                         reset_index(drop=True))
print('MITEs in genes with exons and cds checked:', '\n', mites_in_genes_ex_cds.to_string(max_rows=100))

# duplicated = mites_in_genes_ex_cds[mites_in_genes_ex_cds.duplicated(['chr', 'start_gene', 'end_gene', 'gene_strand',
#                                                                       'unique_no', 'family', 'start_te',
#                                                                       'start_ex', 'end_ex'], keep=False)]
# print('Duplicated:', '\n', duplicated.to_string())

mites_on_board_gene = (bins_to_compare.merge(mites_with_unique_no, how='outer', on=['chr', 'unique_no']).
                       query('(start_te < start_gene) & (end_te >= start_gene) | '
                             '(end_te > end_gene) & (start_te <= end_gene)').
                       reset_index(drop=True))
print('MITEs on genes board:', '\n', mites_on_board_gene.to_string(max_rows=30))

mites_ob_gene_in_exons = (mites_on_board_gene.merge(exons_el10, how='outer', on=['chr']).
                          query('(start_te < start_ex & end_te >= start_ex) | (start_te < end_ex & end_te >= end_ex)').
                          drop_duplicates().reset_index(drop=True))
print('MITEs on genes board in exons:', '\n', mites_ob_gene_in_exons.to_string(max_rows=50))

mites_ob_in_exons = (mites_in_genes.merge(exons_el10, how='outer', on=['chr']).
                     query('(start_te < start_ex & end_te >= start_ex) | (start_te < end_ex & end_te >= end_ex)').
                     drop_duplicates().reset_index(drop=True))
print('MITEs on board in exons:', '\n', mites_ob_in_exons.to_string(max_rows=30))

mites_on_board_genes_exons = (mites_ob_gene_in_exons.merge(mites_ob_in_exons, how='outer').
                              sort_values(by=['chr', 'start_gene', 'start_te']).reset_index(drop=True))
print('MITEs on board in genes and exons:', '\n', mites_on_board_genes_exons.to_string(max_rows=30))

mites_on_board_cds = (mites_on_board_genes_exons.merge(cds_el10_drop, how='outer', on=['chr']).
                      query('(start_te < start_cds & end_te >= start_cds) | (start_te < end_cds & end_te >= end_cds)').
                      drop_duplicates().reset_index(drop=True))
print('MITEs on board in cds:', '\n', mites_on_board_cds.to_string(max_rows=30))

mites_on_board_genes_exons_cds = mites_on_board_genes_exons.merge(mites_on_board_cds, how='outer').reset_index(drop=True)
print('MITEs on board in genes, exons and cds:', '\n', mites_on_board_genes_exons_cds.to_string(max_rows=50))

mites_ob_with_cds_check = ((mites_on_board_genes_exons_cds.merge(cds_el10, how='outer',
                                                               on=['chr', 'gene_name', 'gene_strand']).
                           dropna(subset=['unique_no'])).
                           rename(columns={'start_cds_x': 'start_cds', 'end_cds_x': 'end_cds',
                                           'start_cds_y': 'start_annot_cds', 'end_cds_y': 'end_annot_cds'}))
print('MITEs on board with cds coordinants:', '\n', mites_ob_with_cds_check.to_string(max_rows=50))

x = mites_ob_with_cds_check['start_ex']
y = mites_ob_with_cds_check['start_te']
z = mites_ob_with_cds_check['gene_strand']
g = mites_ob_with_cds_check['start_gene']
c = mites_ob_with_cds_check['start_annot_cds']
e = mites_ob_with_cds_check['end_annot_cds']

# cond_board = [pd.isnull(x) & ((y.lt(g) & z.eq('+')) | (y.gt(g) & z.eq('-'))),
#               pd.isnull(x) & ((y.lt(g) & z.eq('-')) | (y.gt(g) & z.eq('+'))),
#               pd.notnull(x) & ((y.lt(x) & z.eq('+')) | (y.gt(x) & z.eq('-'))),
#               pd.notnull(x) & ((y.lt(x) & z.eq('-')) | (y.gt(x) & z.eq('+')))]
# choices = ['board_up', 'board_down', 'board_exon_up', 'board_exon_down']
# mites_on_board_exons['mite_loc'] = np.select(cond_board, choices, default=0)

cond_utr = [((z.eq('+') & y.lt(c)) | (z.eq('-') & y.gt(c))),
            ((z.eq('-') & y.lt(e)) | (z.eq('+') & y.gt(e)))]

choices_utr = ['5_UTR', '3_UTR']
# mites_on_board_exons['utr_check'] = np.select(cond_utr, choices_utr, default=0)
mites_ob_with_cds_check['mite_loc'] = np.select(cond_utr, choices_utr, default='exon')

mite_loc = mites_ob_with_cds_check.pop('mite_loc')
mites_ob_with_cds_check.insert(11, mite_loc.name, mite_loc)
mites_ob_with_cds_check = mites_ob_with_cds_check.drop_duplicates(
    subset=['chr', 'start_gene', 'end_gene', 'gene_strand', 'unique_no', 'family', 'start_te', 'end_te', 'start_ex',
            'end_ex']).reset_index(drop=True)
print('MITEs on board in genes, exons and checked:', '\n', mites_ob_with_cds_check.to_string())

mites_updown = (bins_to_compare.merge(mites_with_unique_no, how='outer', on=['chr', 'unique_no']).
                query('(start_te >= start_gene - 2000) & (end_te < start_gene) | (end_te <= end_gene + 2000) & '
                      '(start_te > end_gene)').
                reset_index(drop=True))
print('MITEs updown:', '\n', mites_updown.to_string(max_rows=30))

st = mites_updown['start_te']
sg = mites_updown['start_gene']
strand = mites_updown['gene_strand']

cond_updown = [((st.lt(sg) & strand.eq('+')) | (st.gt(sg) & strand.eq('-'))),
               ((st.lt(sg) & strand.eq('-')) | (st.gt(sg) & strand.eq('+')))]

choices_ud = ['upstream', 'downstream']
mites_updown['mite_loc'] = np.select(cond_updown, choices_ud, default=0)
print(mites_updown.to_string(max_rows=30))

exons_board_updown = (pd.concat([mites_in_genes_ex_cds, mites_ob_with_cds_check, mites_updown]).
                      sort_values(by=['chr', 'start_gene', 'start_te', 'mite_loc']).
                      drop_duplicates(subset=['chr', 'start_gene', 'end_gene', 'gene_strand', 'unique_no', 'family',
                      'start_te', 'end_te'], keep='first').reset_index(drop=True))

print('MITEs with <mite_loc>:', '\n', exons_board_updown.to_string(max_rows=400))
cds_check = exons_board_updown[exons_board_updown['mite_loc'].eq('cds')]
print('MITEs cds:', '\n', cds_check)

bins_genes_mites = (opt_bins_and_de.merge(exons_board_updown, how='outer',
                                          on=['chr', 'start_gene', 'end_gene', 'gene_strand', 'unique_no', 'set_no']).
                    reset_index(drop=True))
print(bins_genes_mites.to_string(max_rows=200))

te_ref = pd.read_csv('../files/P4_ref_matrix_sort_2000_all.csv', sep='\t')
te_ref = te_ref.drop(columns=te_ref.columns[-12:], axis=1)
print(te_ref.to_string(max_rows=30))
te_nonref = pd.read_csv('../files/P4_nonref_matrix_all.csv', sep='\t')
te_nonref = te_nonref.drop(columns=te_nonref.columns[-12:], axis=1)
print(te_nonref.to_string(max_rows=30))

ref_check = bins_genes_mites.merge(te_ref, on=['chr', 'family', 'start_te', 'end_te'])
ref_check['te_type'] = 'ref'
print('Ref check:', '\n', ref_check.to_string(max_rows=30))

nonref_check = bins_genes_mites.merge(te_nonref, on=['chr', 'family', 'start_te', 'end_te'])
nonref_check['te_type'] = 'nonref'
print('Nonref check:', '\n', nonref_check.to_string(max_rows=30))

type_checked = pd.concat([ref_check, nonref_check]).reset_index(drop=True)
print('Type checked:', '\n', type_checked.to_string(max_rows=30))

mites_in_genes_te_types = bins_genes_mites.merge(type_checked, how='outer')

te_type = mites_in_genes_te_types.pop('te_type')
mites_in_genes_te_types.insert(45, te_type.name, te_type)
print(mites_in_genes_te_types.to_string(max_rows=300))

# mites_in_genes_te_types.to_csv('../files/P4_all_results_in_one_utr.csv', sep='\t', index=False)
