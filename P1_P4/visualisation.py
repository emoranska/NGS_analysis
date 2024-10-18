import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns

superfam_de_005 = (pd.read_csv('../files/heatmaps/P4_DE_005_superfam_counts.csv', sep='\t').
                   set_index('superfamily'))

superfam_all_genes = (pd.read_csv('../files/heatmaps/P4_all_genes_superfam_counts.csv', sep='\t').
                      set_index('superfamily'))

superfam_all_genes_utr = (pd.read_csv('../files/heatmaps/P4_superfam_all_genes_utr.csv', sep='\t').
                          set_index('superfamily'))

superfam_de_005_utr = (pd.read_csv('../files/heatmaps/P4_superfam_DE005_utr.csv', sep='\t').
                       set_index('superfamily'))

families_de_005 = pd.read_csv('../files/heatmaps/P4_families_DE_005.csv', sep='\t').set_index('family')
families_all_genes = pd.read_csv('../files/heatmaps/P4_families_all_genes.csv', sep='\t').set_index('family')
fam_all_genes_utr = (pd.read_csv('../files/heatmaps/P4_fam_all_genes_utr.csv', sep='\t').set_index('family'))
fam_de_005_utr = pd.read_csv('../files/heatmaps/P4_fam_DE005_utr.csv', sep='\t').set_index('family')

superfam_de_005_updown = (pd.read_csv('../files/heatmaps/P4_DE_005_superfam_updown.csv', sep='\t').
                          set_index('superfamily'))

superfam_updown_utr = (pd.read_csv('../files/heatmaps/P4_superfam_updown_utr.csv', sep='\t').
                       set_index('superfamily'))

families_de_005_updown = (pd.read_csv('../files/heatmaps/P4_DE_005_families_updown.csv', sep='\t').
                          set_index('family'))

fam_updown_utr = (pd.read_csv('../files/heatmaps/P4_fam_updown_utr.csv', sep='\t').set_index('family'))


def heatmap_mites_loc(data, title, ylabel_fontsize, values_fontsize, cbar_param):
    # simple heatmap in Excel
    # df.style.background_gradient(cmap='Blues').set_properties(**{'font-size': '20px'}).to_excel('styled.xlsx')

    plt.figure(figsize=(10, 100))

    # 'annot' and 'square' params to determine according to input data, cmap=YlOrBr/YlGnBu
    hm = sns.heatmap(data, cmap="YlOrBr", annot=False, fmt="0.0f", annot_kws={'size': values_fontsize}, linewidths=.5,
                     cbar_kws=cbar_param)
    plt.title(title, fontsize=20)
    plt.xlabel('Localisation', fontsize=12)
    plt.ylabel('Family', fontsize=12)  # or 'Superfamily'

    cbar = hm.collections[0].colorbar
    cbar.ax.tick_params(labelsize=9)

    hm.tick_params(axis='y', labelsize=ylabel_fontsize)
    hm.tick_params(axis='x', labelsize=9, rotation=50)

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

# heatmap_mites_loc(superfam_all_genes_utr[['upstream', '5_UTR', 'intron', 'cds', '3_UTR', 'downstream']],
#                   'P4 MITEs (for all genes)',
#        10, 10, {"pad": 0.02, "shrink": 0.5})

# heatmap_mites_loc(fam_all_genes_utr[['upstream', '5_UTR', 'intron', 'cds', '3_UTR', 'downstream']],
#                   'P4 MITEs (for all genes)',5, 5, {"pad": 0.01, "shrink": 0.35})

# heatmap_mites_loc(superfam_de_005_utr[['upstream', '5_UTR', 'intron', '3_UTR', 'downstream']],
#                   'P4 MITEs (for DE genes < 0.05)',10, 10, {"pad": 0.02, "shrink": 0.5})

# heatmap_mites_loc(fam_de_005_utr[['upstream', '5_UTR', 'intron', 'cds', '3_UTR', 'downstream']],
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

# to check and change!
def heatmap_subplots(data, title, ylabel_fontsize, values_fontsize, cbar_param):
    # fig, axes = plt.subplots(1, 2, sharex='all', figsize=(10, 5))
    # fig.suptitle(title)
    # axes[0].set_title('up-regulated')
    # axes[1].set_title('down-regulated')

    up_reg = data[data.columns[:5]]
    down_reg = data[data.columns[5:]]

    fig, (ax1, ax2) = plt.subplots(ncols=2, figsize=(10, 100))
    fig.suptitle(title)

    ax1.set_title('up-regulated')
    # ax2.set_title('down-regulated')

    hm1 = sns.heatmap(up_reg, ax=ax1, cmap="YlOrBr", cbar=False, linewidths=.1, )
    hm2 = sns.heatmap(down_reg, ax=ax2, cmap="YlOrBr", cbar_kws=cbar_param, linewidths=.1, )
    hm2.set(yticklabels=[])
    hm2.set(title='down-regulated')
    hm2.set(ylabel=None)

    hm1.tick_params(axis='y', labelsize=ylabel_fontsize)
    hm1.tick_params(axis='x', labelsize=8, rotation=45)
    hm2.tick_params(axis='x', labelsize=8, rotation=45)
    hm2.tick_params(left=False)

    # plt.figure(figsize=(10, 100))
    # hm = sns.heatmap(data, cmap="YlOrBr", annot=True, fmt="0.0f", annot_kws={'size': values_fontsize}, linewidths=.5,
    #                  cbar_kws=cbar_param, square=True)
    # plt.title(title, fontsize=20)
    # plt.xlabel('Localisation', fontsize=12)
    # plt.ylabel('Superfamily', fontsize=12)
    #
    # cbar = hm.collections[0].colorbar
    # cbar.ax.tick_params(labelsize=8)
    #
    # hm.tick_params(axis='y', labelsize=ylabel_fontsize)
    plt.show()


# heatmap_subplots(superfam_de_005_updown[superfam_de_005_updown.columns[2:]], 'P4 MITEs associated with DEGs',
#                9, 7, {"pad": 0.01, "shrink": 0.35})

# no cds
heatmap_subplots(superfam_updown_utr[['upstream_up', '5_UTR_up', 'intron_up', '3_UTR_up', 'downstream_up','upstream_down',
                                  '5_UTR_down', 'intron_down', '3_UTR_down', 'downstream_down']],
                  'P4 MITEs associated with DEGs',
                  9, 5, {"pad": 0.01, "shrink": 0.35})

# heatmap_subplots(fam_updown_utr[['upstream_up', '5_UTR_up', 'intron_up', '3_UTR_up', 'downstream_up','upstream_down',
#                                   '5_UTR_down', 'intron_down', '3_UTR_down', 'downstream_down']],
#                   'P4 MITEs associated with DEGs',
#                   6, 5, {"pad": 0.01, "shrink": 0.35})

# with cds
# heatmap_subplots(fam_updown_utr[['upstream_up', '5_UTR_up', 'intron_up', 'cds_up', '3_UTR_up', 'downstream_up','upstream_down',
#                                   '5_UTR_down', 'intron_down', 'cds_down', '3_UTR_down', 'downstream_down']],
#                   'P4 MITEs associated with DEGs',
#                   6, 5, {"pad": 0.01, "shrink": 0.35})
