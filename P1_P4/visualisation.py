import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns

superfam_de_005 = (pd.read_csv('../files/P4_DE_005_superfam_counts.csv', sep='\t').
                   set_index('superfamily'))

superfam_all_genes = (pd.read_csv('../files/P4_all_genes_superfam_counts.csv', sep='\t').
                      set_index('superfamily'))

families_de_005 = pd.read_csv('../files/P4_families_DE_005.csv', sep='\t').set_index('family')
families_all_genes = pd.read_csv('../files/P4_families_all_genes.csv', sep='\t').set_index('family')

superfam_de_005_updown = (pd.read_csv('../files/P4_DE_005_superfam_updown.csv', sep='\t').
                          set_index('superfamily'))

families_de_005_updown = (pd.read_csv('../files/P4_DE_005_families_updown.csv', sep='\t').
                          set_index('family'))


def heatmap_mites_loc(data, title, ylabel_fontsize, values_fontsize, cbar_param):
    # simple heatmap in Excel
    # df.style.background_gradient(cmap='Blues').set_properties(**{'font-size': '20px'}).to_excel('styled.xlsx')

    plt.figure(figsize=(10, 100))

    # 'annot' and 'square' params to determine according to input data
    hm = sns.heatmap(data, cmap="YlOrBr", annot=False, fmt="0.0f", annot_kws={'size': values_fontsize}, linewidths=.5,
                     cbar_kws=cbar_param)
    plt.title(title, fontsize=20)
    plt.xlabel('Localisation', fontsize=12)
    plt.ylabel('Family', fontsize=12)  # or 'Family'

    cbar = hm.collections[0].colorbar
    cbar.ax.tick_params(labelsize=9)

    hm.tick_params(axis='y', labelsize=ylabel_fontsize)
    hm.tick_params(axis='x', labelsize=9, rotation=45)

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


# to check and change!
def heatmap_subplots(data, title, ylabel_fontsize, values_fontsize, cbar_param):
    fig, axes = plt.subplots(1, 2, sharex='all', figsize=(10, 5))
    fig.suptitle(title)
    axes[0].set_title('up-regulated')
    axes[1].set_title('down-regulated')

    up_reg = data[data.columns[2:5]]
    down_reg = data[data[6:9]]

    fig, (ax1, ax2) = plt.subplots(ncols=2)

    sns.heatmap(up_reg, ax=ax1)
    sns.heatmap(down_reg, ax=ax2)

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


heatmap_subplots(superfam_de_005_updown[superfam_de_005_updown.columns[2:]], 'P4 MITEs associated with DEGs',
               9, 7, {"pad": 0.01, "shrink": 0.35})
