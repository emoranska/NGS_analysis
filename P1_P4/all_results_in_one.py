import pandas as pd
import numpy as np
# from concat_bins_with_genes_and_rest import find_unique_no

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

del bins_and_de
del bins_and_de_obj
del bins_conv_obj
del bins_and_de_float
del bins_conv_float
del bins_and_de_int
del bins_conv_int

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

del mites_all_obj, mites_conv_obj, mites_all_float, mites_conv_float, mites_all_int, mites_conv_int

mites_with_unique_no = (((opt_mites_all.merge(sets_list, how='cross').
                        query('((te_ins == one_one) & (no_te_ins == zero_zero)) | '
                              '((te_ins == zero_zero) & (no_te_ins == one_one))')).
                        drop(columns=['one_one', 'zero_zero']).rename(columns={'start': 'start_te', 'end': 'end_te'})).
                        reset_index(drop=True))

# opt_mites_all = None
del opt_mites_all
print('MITEs with unique_no:', '\n', mites_with_unique_no.to_string(max_rows=30))

# check MITEs in genes and on board's genes
mites_in_genes = (bins_to_compare.merge(mites_with_unique_no, how='outer', on=['chr', 'unique_no']).
                  query('((start_te >= start_gene) & (end_te <= end_gene)) | '
                        '((start_te < start_gene) & (end_te >= start_gene)) | '
                        '((end_te > end_gene) & (start_te <= end_gene))').drop_duplicates().reset_index(drop=True))
mites_in_genes['mite_loc'] = 'intron'
print('MITEs in genes and on board genes:', '\n', mites_in_genes.to_string(max_rows=30))

# check MITEs in exons and on board's exons
exons_el10 = pd.read_csv('../files/EL10_exons_longest_isoform.csv', sep='\t')

opt_exons_el10 = exons_el10.copy()
exons_el10_obj = exons_el10.select_dtypes(include=['object']).copy()

exons_el10_conv_obj = pd.DataFrame()
for col in exons_el10_obj.columns:
    num_unique_values = len(exons_el10_obj[col].unique())
    num_total_values = len(exons_el10_obj[col])
    if num_unique_values / num_total_values < 0.5:
        exons_el10_conv_obj.loc[:, col] = exons_el10_obj[col].astype('category')
    else:
        exons_el10_conv_obj.loc[:, col] = exons_el10_obj[col]

opt_exons_el10[exons_el10_conv_obj.columns] = exons_el10_conv_obj
exons_el10_int = exons_el10.select_dtypes(include=['int'])
exons_el10_conv_int = exons_el10_int.apply(pd.to_numeric, downcast='unsigned')
opt_exons_el10[exons_el10_conv_int.columns] = exons_el10_conv_int

del exons_el10, exons_el10_obj, exons_el10_conv_obj, exons_el10_int, exons_el10_conv_int

mites_in_exons = (mites_in_genes.merge(opt_exons_el10, how='outer', on=['chr', 'gene_name', 'gene_strand']).
                  query('(start_te >= start_ex & end_te <= end_ex) | '
                        '(start_te < start_ex & end_te >= start_ex) | '
                        '(start_te < end_ex & end_te >= end_ex)').
                  drop_duplicates(subset=['chr', 'start_gene', 'end_gene', 'gene_name', 'gene_strand', 'unique_no',
                                          'family', 'start_te', 'end_te']).reset_index(drop=True))
mites_in_exons['mite_loc'] = 'exon'
print('MITEs in exons and on board exons:', '\n', mites_in_exons.to_string(max_rows=30))

# find MITEs in cds
cds_el10 = pd.read_csv('../files/EL10_CDS_bed.csv', sep='\t')

mites_in_genes_with_cds_check = (mites_in_genes.merge(cds_el10, how='outer', on=['chr', 'gene_name', 'gene_strand']).
                                 dropna(subset=['unique_no']).query('(start_te >= start_cds & end_te <= end_cds) | '
                                                                    '(start_te < start_cds & end_te >= start_cds) | '
                                                                    '(start_te <= end_cds & end_te > end_cds)').
                                 drop_duplicates().reset_index(drop=True))

mites_in_genes_with_cds_check['mite_loc'] = 'cds'
print('MITEs in cds:', '\n', mites_in_genes_with_cds_check.to_string(max_rows=50))

# find MITEs in UTRs
el10_utrs = pd.read_csv('../files/EL10_UTRs_longest_isoform.csv', sep='\t')
mites_in_genes_with_utr_check = (mites_in_genes.merge(el10_utrs, on=['chr', 'gene_name', 'gene_strand']).
                                 query('(start_te >= start_utr & end_te <= end_utr) | '
                                       '(start_te < start_utr & end_te > start_utr) | '
                                       '(start_te < end_utr & end_te > end_utr)').drop_duplicates().
                                 reset_index(drop=True)).drop(columns=['mite_loc']).rename(columns={'utr': 'mite_loc'})
print('MITEs in UTRs:', '\n', mites_in_genes_with_utr_check.to_string())

# concat MITEs in cds and UTRs
cds_and_utrs = ((pd.concat([mites_in_genes_with_cds_check, mites_in_genes_with_utr_check]).
                sort_values(by=['chr', 'start_gene', 'start_te', 'unique_no', 'mite_loc'])).
                drop_duplicates(subset=['chr', 'start_gene', 'end_gene', 'gene_name', 'gene_strand', 'unique_no',
                                        'family', 'start_te', 'end_te'], keep='last').reset_index(drop=True))
print('MITEs in cds and UTRs together', '\n', cds_and_utrs.to_string())

# MITEs in cds, UTRs and introns merged together
mites_in_genes_with_cds_and_utrs = (mites_in_genes.merge(cds_and_utrs, how='outer').
                                    sort_values(by=['chr', 'start_gene', 'start_te', 'mite_loc']).
                                    drop_duplicates(subset=['chr', 'start_gene', 'end_gene', 'gene_name', 'gene_strand',
                                                            'unique_no', 'family', 'start_te', 'end_te']).
                                    reset_index(drop=True))
print('MITEs in genes with cds and UTRs together', '\n', mites_in_genes_with_cds_and_utrs.to_string(max_rows=100))

cds_check = mites_in_genes_with_cds_and_utrs[mites_in_genes_with_cds_and_utrs['mite_loc'].eq('cds')]
print('MITEs cds:', '\n', cds_check)

utr5_check = mites_in_genes_with_cds_and_utrs[mites_in_genes_with_cds_and_utrs['mite_loc'].eq('5UTR')]
print('MITEs 5UTR:', '\n', utr5_check)

# check exons not in cds with lncRNA reference annotation
# exon_no_cds_check = mites_in_genes_with_cds_check[mites_in_genes_with_cds_check['mite_loc'].eq('exon')]
# # print('MITEs in exons not in cds:', '\n', exon_no_cds_check.reset_index(drop=True).to_string())
#
# lncRNA_el10 = pd.read_csv('../files/EL10_lncRNA.csv', sep='\t')
# lncRNA_check = exon_no_cds_check.merge(lncRNA_el10, how='inner', on=['chr', 'start_ex', 'end_ex'])
# print('MITEs in lncRNA regions:', '\n', lncRNA_check.to_string())
# mites_in_genes_ex.loc[~np.isnan(mites_in_genes_ex['start_ex']), 'mite_loc'] = 'exon'

# find MITEs upstream and downstream +/- 2000 bp
mites_updown = (bins_to_compare.merge(mites_with_unique_no, how='outer', on=['chr', 'unique_no']).
                query('((start_te >= start_gene - 2000) & (end_te < start_gene)) | ((end_te <= end_gene + 2000) & '
                      '(start_te > end_gene))').
                reset_index(drop=True))
print('MITEs updown:', '\n', mites_updown.to_string(max_rows=30))

del mites_with_unique_no

st = mites_updown['start_te']
et = mites_updown['end_te']
sg = mites_updown['start_gene']
eg = mites_updown['end_gene']
strand = mites_updown['gene_strand']

cond_updown = [((et.lt(sg) & strand.eq('+')) | (st.gt(sg) & strand.eq('-'))),
               ((st.lt(sg) & strand.eq('-')) | (st.gt(eg) & strand.eq('+')))]

choices_ud = ['upstream', 'downstream']
mites_updown['mite_loc'] = np.select(cond_updown, choices_ud, default=0)
print(mites_updown.to_string(max_rows=30))

upstream_check = mites_updown[mites_updown['mite_loc'].eq('upstream')]
print('MITEs upstream check:', '\n', upstream_check.reset_index(drop=True).to_string(max_rows=50))

downstream_check = mites_updown[mites_updown['mite_loc'].eq('downstream')]
print('MITEs downstream check:', '\n', downstream_check.reset_index(drop=True).to_string(max_rows=50))

exons_board_updown = (pd.concat([mites_in_genes_with_cds_and_utrs, mites_updown]).
                      sort_values(by=['chr', 'start_gene', 'start_te', 'mite_loc']).
                      drop_duplicates(subset=['chr', 'start_gene', 'end_gene', 'gene_strand', 'unique_no', 'family',
                      'start_te', 'end_te', 'mite_loc'], keep='first').reset_index(drop=True))

print('MITEs with <mite_loc>:', '\n', exons_board_updown.to_string(max_rows=200))
cds_check = exons_board_updown[exons_board_updown['mite_loc'].eq('cds')]
print('MITEs cds:', '\n', cds_check)

bins_genes_mites = (opt_bins_and_de.merge(exons_board_updown, how='outer',
                                          on=['chr', 'gene_name', 'start_gene', 'end_gene', 'gene_strand', 'unique_no', 'set_no']).
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

mites_in_genes_te_types.to_csv('../files/P4_all_results_in_one_utr_longest_isoforms.csv', sep='\t', index=False)
