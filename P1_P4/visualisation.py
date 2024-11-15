import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns

superfam_de_005 = (pd.read_csv('../files/heatmaps/hm_data/P4_DE_005_superfam_counts.csv', sep='\t').
                   set_index('superfamily'))

superfam_all_genes = (pd.read_csv('../files/heatmaps/hm_data/P4_all_genes_superfam_counts.csv', sep='\t').
                      set_index('superfamily'))

superfam_all_genes_utr = (pd.read_csv('../files/heatmaps/hm_data/P4_superfam_all_genes_st.csv', sep='\t',
                                      dtype={'upstream': np.float64, '5_UTR': np.float64, 'intron': np.float64,
                                             'cds': np.float64, '3_UTR': np.float64, 'downstream': np.float64}).
                          set_index('superfamily'))

superfam_de_005_utr = (pd.read_csv('../files/heatmaps/hm_data/P4_superfam_DE005_st.csv', sep='\t').
                       set_index('superfamily'))

families_de_005 = (pd.read_csv('../files/heatmaps/hm_data/P4_families_DE_005.csv', sep='\t').
                   set_index('family'))
families_all_genes = (pd.read_csv('../files/heatmaps/hm_data/P4_families_all_genes.csv', sep='\t').
                      set_index('family'))
fam_all_genes_utr = (pd.read_csv('../files/heatmaps/hm_data/P4_fam_all_genes_st.csv', sep='\t').
                     set_index('family'))
fam_de_005_utr = (pd.read_csv('../files/heatmaps/hm_data/P4_fam_DE005_st.csv', sep='\t').
                  set_index('family'))

superfam_de_005_updown = (pd.read_csv('../files/heatmaps/hm_data/P4_DE_005_superfam_updown.csv', sep='\t').
                          set_index('superfamily'))

superfam_updown_utr = (pd.read_csv('../files/heatmaps/hm_data/P4_superfam_updown_st.csv', sep='\t').
                       set_index('superfamily'))

families_de_005_updown = (pd.read_csv('../files/heatmaps/hm_data/P4_DE_005_families_updown.csv', sep='\t').
                          set_index('family'))

fam_updown_utr = (pd.read_csv('../files/heatmaps/hm_data/P4_fam_updown_st.csv', sep='\t').
                  set_index('family'))

fam_6 = pd.read_csv('../files/heatmaps/hm_data/6_families_DEGs.csv', sep='\t').set_index('family')
fam_6_st = pd.read_csv('../files/heatmaps/hm_data/6_families_DEGs_st.csv', sep='\t').set_index('family')

fam_6_updown = pd.read_csv('../files/heatmaps/hm_data/6_fam_DEGs_updown.csv', sep='\t').set_index('family')
fam_6_updown_st = (pd.read_csv('../files/heatmaps/hm_data/6_fam_DEGs_updown_st.csv', sep='\t').
                   set_index('family'))

fam_6_all_st = pd.read_csv('../files/heatmaps/hm_data/6_families_all_st.csv', sep='\t').set_index('family')


def heatmap_mites_loc(data, title, ylabel_fontsize, values_fontsize, cbar_param):
    # simple heatmap in Excel
    # df.style.background_gradient(cmap='Blues').set_properties(**{'font-size': '20px'}).to_excel('styled.xlsx')

    plt.figure(figsize=(9, 10))

    # 'annot' and 'square' params to determine according to input data, cmap=YlOrBr/YlGnBu
    hm = sns.heatmap(data, cmap="YlOrBr", annot=False, fmt="0.0f", annot_kws={'size': values_fontsize}, linewidths=.5,
                     cbar_kws=cbar_param, square=True)
    plt.title(title, fontsize=20)
    plt.xlabel('Localisation', fontsize=12)
    plt.ylabel('Family', fontsize=12)  # or 'Superfamily'

    cbar = hm.collections[0].colorbar
    cbar.ax.tick_params(labelsize=9)

    hm.tick_params(axis='y', labelsize=ylabel_fontsize)
    hm.tick_params(axis='x', labelsize=9, rotation=50)

    hm_figure = hm.get_figure()
    hm_figure.savefig(f'../files/heatmaps/hm_6_families/6_families.pdf')
    plt.show()


# heatmap_mites_loc(superfam_de_005[['upstream', 'downstream', 'intron', 'cds']],
# 'P4 MITEs (for DE genes < 0.05)', 10, 10, {"pad": 0.02, "shrink": 0.5})
# heatmap_mites_loc(superfam_all_genes[['upstream', 'downstream', 'intron', 'cds']], 'P4 MITEs (for all genes)',
#        10, 10, {"pad": 0.02, "shrink": 0.5})
# heatmap_mites_loc(families_de_005[['upstream', 'downstream', 'intron', 'cds']], 'P4 MITEs (for DE genes < 0.05)',
#                   6, 7, {"pad": 0.01, "shrink": 0.35})
# heatmap_mites_loc(families_all_genes[['upstream', 'downstream', 'intron', 'cds']], 'P4 MITEs (for all genes)',
#         5,5, {"pad": 0.01, "shrink": 0.35})

# heatmap_mites_loc(superfam_de_005_updown[superfam_de_005_updown.columns[2:]], 'P4 MITEs associated with DEGs',
#                   9, 7, {"pad": 0.01, "shrink": 0.35})

# heatmap_mites_loc(families_de_005_updown, 'P4 MITEs associated with DEGs', 6, 5, {"pad": 0.01, "shrink": 0.35})

# heatmap_mites_loc(superfam_all_genes_utr,
#                   'P4 MITEs (for all genes)',
#        10, 10, {"pad": 0.02, "shrink": 0.5})

# heatmap_mites_loc(fam_all_genes_utr[['upstream', '5_UTR', 'intron', '3_UTR', 'downstream']],
#                   'P4 MITEs (for all genes)',5, 5, {"pad": 0.01, "shrink": 0.35})

# heatmap_mites_loc(superfam_de_005_utr,
#                   'P4 MITEs (for DE genes < 0.05)',10, 10, {"pad": 0.02, "shrink": 0.5})

# heatmap_mites_loc(fam_de_005_utr[['upstream', '5_UTR', 'intron', '3_UTR', 'downstream']],
#                   'P4 MITEs (for DE genes < 0.05)', 6, 5, {"pad": 0.01, "shrink": 0.35})

# heatmap_mites_loc(superfam_updown_utr[superfam_updown_utr.columns[3:]], 'P4 MITEs associated with DEGs',
#                   9, 7, {"pad": 0.01, "shrink": 0.35})

# heatmap_mites_loc(superfam_updown_utr[['upstream_up', '5_UTR_up', 'intron_up', '3_UTR_up', 'downstream_up',
#                                        'upstream_down', '5_UTR_down', 'intron_down', '3_UTR_down', 'downstream_down']],
#                   'P4 MITEs associated with DEGs',
#                   9, 7, {"pad": 0.01, "shrink": 0.35})

# heatmap_mites_loc(fam_updown_utr[fam_updown_utr.columns[3:]], 'P4 MITEs associated with DEGs',
#                   6, 5, {"pad": 0.01, "shrink": 0.35})

# heatmap_mites_loc(fam_updown_utr[['upstream_up', '5_UTR_up', 'intron_up', '3_UTR_up', 'downstream_up','upstream_down',
#                                   '5_UTR_down', 'intron_down', '3_UTR_down', 'downstream_down']],
#                   'P4 MITEs associated with DEGs',
#                   6, 5, {"pad": 0.01, "shrink": 0.35})

# heatmap_mites_loc(fam_6[['upstream', '5_UTR', 'intron', '3_UTR', 'downstream']], '6 MITEs families',
#                   10, 10, {"pad": 0.02, "shrink": 0.5})

# heatmap_mites_loc(fam_6_st[['upstream', '5_UTR', 'intron', '3_UTR', 'downstream']],
#                   '6 MITEs families - standarised (per 10 Mb)',
#                   10, 10, {"pad": 0.02, "shrink": 0.5})

heatmap_mites_loc(fam_6_all_st[['upstream', '5_UTR', 'intron', '3_UTR', 'downstream']],
                  '6 MITEs families (all genes) - standarised (per 10 Mb)',
                  10, 10, {"pad": 0.02, "shrink": 0.5})


def heatmap_subplots(data, title, ylabel_fontsize):

    # with cds
    # up_reg = data[data.columns[:6]]
    # down_reg = data[data.columns[6:]]

    # without cds
    up_reg = data[data.columns[:5]]
    down_reg = data[data.columns[5:]]

    # for superfamilies
    # fig, (ax1, ax2) = plt.subplots(ncols=2, figsize=(12, 6), sharey='all')

    # for families
    fig, (ax1, ax2) = plt.subplots(ncols=2, figsize=(10, 100), sharey='all')

    fig.suptitle(title)

    # labels_x = ['upstream', '5_UTR', 'intron', 'cds', '3_UTR', 'downstream']

    # without cds
    labels_x = ['upstream', '5_UTR', 'intron', '3_UTR', 'downstream']

# vmax = 4.0 for families, 6.0 for superfamilies, 8.0 for 6_fam_updown, 2.0 for 6_fam_updown_st
    hm1 = sns.heatmap(up_reg, ax=ax1, cmap="YlOrBr", cbar=False, linewidths=.1, xticklabels=labels_x, vmin=0, vmax=2.0)
    hm2 = sns.heatmap(down_reg, ax=ax2, cmap="YlOrBr", cbar=False, linewidths=.1, xticklabels=labels_x, vmin=0, vmax=2.0)

    hm1.set(title='up-regulated')
    hm2.set(title='down-regulated')
    hm2.set(ylabel=None)

    hm1.tick_params(axis='y', labelsize=ylabel_fontsize)
    hm1.tick_params(axis='x', labelsize=8, rotation=45)
    hm2.tick_params(axis='x', labelsize=8, rotation=45)
    hm2.tick_params(left=False)

    im = plt.gca().get_children()[0]
    # cax = fig.add_axes([0.93, 0.3, 0.01, 0.35])  # for families
    cax = fig.add_axes([0.93, 0.3, 0.01, 0.35])  # for superfamilies
    fig.colorbar(im, cax=cax)

    plt.show()


# with cds
# heatmap_subplots(superfam_updown_utr, '4 MITEs associated with DEGs', 9)

# heatmap_subplots(fam_updown_utr, 'P4 MITEs associated with DEGs', 6)

# no cds
# heatmap_subplots(fam_updown_utr[['upstream_up', '5_UTR_up', 'intron_up', '3_UTR_up', 'downstream_up', 'upstream_down',
#                                  '5_UTR_down', 'intron_down', '3_UTR_down', 'downstream_down']],
#                  'P4 MITEs associated with DEGs', 6,)

# heatmap_subplots(fam_6_updown[['upstream_up', '5_UTR_up', 'intron_up', '3_UTR_up', 'downstream_up', 'upstream_down',
#                                  '5_UTR_down', 'intron_down', '3_UTR_down', 'downstream_down']],
#                  '6 MITEs families associated with DEGs', 9)

# heatmap_subplots(fam_6_updown_st[['upstream_up', '5_UTR_up', 'intron_up', '3_UTR_up', 'downstream_up', 'upstream_down',
#                                  '5_UTR_down', 'intron_down', '3_UTR_down', 'downstream_down']],
#                  '6 MITEs families associated with DEGs - standarised (per 10 Mb)', 9)
