import pandas as pd
import numpy as np

# mites_in_el10 = pd.read_csv('../files/MITEs_in_EL10_from_bed.csv', sep='\t')
mites_in_el10 = pd.read_csv('../files/P4_all_MITEs_from_matrix.csv', sep='\t')
genes_el10 = pd.read_csv('../files/EL10_genes_annotation.csv', sep='\t')
cds_el10 = pd.read_csv('../files/EL10_CDS_bed.csv', sep='\t')
utr_el10 = pd.read_csv('../files/EL10_UTRs_longest_isoform.csv', sep='\t')

print(mites_in_el10.memory_usage(deep=True))

mites_in_el10_int = mites_in_el10.select_dtypes(include=['int'])
print(mites_in_el10_int.memory_usage(deep=True))
mites_in_el10_conv_int = mites_in_el10_int.apply(pd.to_numeric, downcast='unsigned')

mites_in_el10_obj = mites_in_el10.select_dtypes(include=['object'])
print(mites_in_el10_obj.memory_usage(deep=True))

mites_in_el10_conv_obj = pd.DataFrame()
for col in mites_in_el10_obj.columns:
    num_unique_values = len(mites_in_el10_obj[col].unique())
    num_total_values = len(mites_in_el10_obj[col])
    if num_unique_values / num_total_values < 0.5:
        mites_in_el10_conv_obj.loc[:, col] = mites_in_el10_obj[col].astype('category')
    else:
        mites_in_el10_conv_obj.loc[:, col] = mites_in_el10_obj[col]

opt_mites_in_el10 = mites_in_el10.copy()
opt_mites_in_el10[mites_in_el10_conv_int.columns] = mites_in_el10_conv_int
opt_mites_in_el10[mites_in_el10_conv_obj.columns] = mites_in_el10_conv_obj
print(opt_mites_in_el10.memory_usage(deep=True))

del mites_in_el10
del mites_in_el10_int
del mites_in_el10_conv_int
del mites_in_el10_obj
del mites_in_el10_conv_obj

print(genes_el10.memory_usage(deep=True))

genes_el10_int = genes_el10.select_dtypes(include=['int'])
print(genes_el10_int.memory_usage(deep=True))
genes_el10_conv_int = genes_el10_int.apply(pd.to_numeric, downcast='unsigned')

genes_el10_obj = genes_el10.select_dtypes(include=['object'])
print(genes_el10_obj.memory_usage(deep=True))

genes_el10_conv_obj = pd.DataFrame()
for col in genes_el10_obj.columns:
    num_unique_values = len(genes_el10_obj[col].unique())
    num_total_values = len(genes_el10_obj[col])
    if num_unique_values / num_total_values < 0.5:
        genes_el10_conv_obj.loc[:, col] = genes_el10_obj[col].astype('category')
    else:
        genes_el10_conv_obj.loc[:, col] = genes_el10_obj[col]

opt_genes_el10 = genes_el10.copy()
opt_genes_el10[genes_el10_conv_int.columns] = genes_el10_conv_int
opt_genes_el10[genes_el10_conv_obj.columns] = genes_el10_conv_obj
print(opt_genes_el10.memory_usage(deep=True))

del genes_el10
del genes_el10_int
del genes_el10_conv_int
del genes_el10_obj
del genes_el10_conv_obj

# check MITEs in genes and on board's genes
mites_in_genes = (opt_genes_el10.merge(opt_mites_in_el10, how='outer', on=['chr']).
                  query('((start_te >= start_gene) & (end_te <= end_gene)) | '
                        '((start_te < start_gene) & (end_te >= start_gene)) | '
                        '((end_te > end_gene) & (start_te <= end_gene))').drop_duplicates().reset_index(drop=True))
mites_in_genes['mite_loc'] = 'intron'
print('MITEs in EL10 genes and on board genes:', '\n', mites_in_genes.to_string(max_rows=30))

# del opt_mites_in_el10

# find MITEs in cds
mites_in_genes_with_cds_check = (mites_in_genes.merge(cds_el10, how='outer', on=['chr', 'gene_name', 'gene_strand']).
                                 query('(start_te >= start_cds & end_te <= end_cds) | '
                                       '(start_te < start_cds & end_te >= start_cds) | '
                                       '(start_te <= end_cds & end_te > end_cds)').
                                 drop_duplicates().reset_index(drop=True))

mites_in_genes_with_cds_check['mite_loc'] = 'cds'
print('MITEs in cds:', '\n', mites_in_genes_with_cds_check.to_string(max_rows=50))

del cds_el10

# find MITEs in UTRs
mites_in_genes_with_utr_check = (mites_in_genes.merge(utr_el10, on=['chr', 'gene_name', 'gene_strand']).
                                 query('(start_te >= start_utr & end_te <= end_utr) | '
                                       '(start_te < start_utr & end_te > start_utr) | '
                                       '(start_te < end_utr & end_te > end_utr)').drop_duplicates().
                                 reset_index(drop=True)).drop(columns=['mite_loc']).rename(columns={'utr': 'mite_loc'})
print('MITEs in UTRs:', '\n', mites_in_genes_with_utr_check.to_string(max_rows=50))

del utr_el10

# concat MITEs in cds and UTRs
cds_and_utrs = ((pd.concat([mites_in_genes_with_cds_check, mites_in_genes_with_utr_check]).
                sort_values(by=['chr', 'start_gene', 'start_te', 'mite_loc'])).
                drop_duplicates(subset=['chr', 'start_gene', 'end_gene', 'gene_name', 'gene_strand',
                                        'family', 'start_te', 'end_te'], keep='last').reset_index(drop=True))
print('MITEs in cds and UTRs together', '\n', cds_and_utrs.to_string(max_rows=50))

cds_check = cds_and_utrs[cds_and_utrs['mite_loc'].eq('cds')]
print('MITEs cds:', '\n', cds_check)

utr5_check = cds_and_utrs[cds_and_utrs['mite_loc'].eq('5UTR')]
print('MITEs 5UTR:', '\n', utr5_check)

utr3_check = cds_and_utrs[cds_and_utrs['mite_loc'].eq('3UTR')]
print('MITEs 3UTR:', '\n', utr3_check)

# MITEs in cds, UTRs and introns merged together
mites_in_genes_with_cds_and_utrs = (mites_in_genes.merge(cds_and_utrs, how='outer').
                                    sort_values(by=['chr', 'start_gene', 'start_te', 'mite_loc']).
                                    drop_duplicates(subset=['chr', 'start_gene', 'end_gene', 'gene_name', 'gene_strand',
                                                            'family', 'start_te', 'end_te']).
                                    reset_index(drop=True))
print('MITEs in genes with cds and UTRs together', '\n', mites_in_genes_with_cds_and_utrs.to_string(max_rows=30))

del cds_and_utrs
del utr5_check, utr3_check, cds_check

# find MITEs upstream and downstream +/- 2000 bp
mites_updown = (opt_genes_el10.merge(opt_mites_in_el10, how='outer', on=['chr']).
                query('((start_te >= start_gene - 2000) & (end_te < start_gene)) | ((end_te <= end_gene + 2000) & '
                      '(start_te > end_gene))').
                reset_index(drop=True))
print('MITEs updown:', '\n', mites_updown.to_string(max_rows=30))

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

del opt_mites_in_el10
del opt_genes_el10

# upstream_check = mites_updown[mites_updown['mite_loc'].eq('upstream')]
# print('MITEs upstream check:', '\n', upstream_check.reset_index(drop=True).to_string(max_rows=50))
#
# downstream_check = mites_updown[mites_updown['mite_loc'].eq('downstream')]
# print('MITEs downstream check:', '\n', downstream_check.reset_index(drop=True).to_string(max_rows=50))

all_mites_in_el10 = (pd.concat([mites_in_genes_with_cds_and_utrs, mites_updown]).
                      sort_values(by=['chr', 'start_gene', 'start_te', 'mite_loc']).
                     drop_duplicates(subset=['chr', 'start_gene', 'end_gene', 'gene_strand', 'family',
                      'start_te', 'end_te', 'mite_loc'], keep='first').
reset_index(drop=True))

all_mites_in_el10.to_csv('../files/P4_genes_with_MITEs.csv', sep='\t', index=False)
print('MITEs with <mite_loc>:', '\n', all_mites_in_el10.to_string(max_rows=200))
