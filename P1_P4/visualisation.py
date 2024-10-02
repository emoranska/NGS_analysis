import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns

superfam_de_005 = (pd.read_csv('../files/P4_DE_005_superfam_counts.csv', sep='\t').
                   set_index('superfamily'))

superfam_all_genes = (pd.read_csv('../files/P4_all_genes_superfam_counts.csv', sep='\t').
                      set_index('superfamily'))

families_de_005 = pd.read_csv('../files/P4_families_DE_005.csv', sep='\t').set_index('family')
families_all_genes = pd.read_csv('../files/P4_families_all_genes.csv', sep='\t').set_index('family')

def heatmap(data, title, ylabel_fontsize, values_fontsize, cbar_param):
    # działa
    # df.style.background_gradient(cmap='Blues').set_properties(**{'font-size': '20px'}).to_excel('styled.xlsx')

    # fig, ax = plt.subplots(figsize=(20, 150))
    plt.figure(figsize=(10, 100))
    hm = sns.heatmap(data, cmap="YlOrBr", annot=True, fmt="0.0f", annot_kws={'size': values_fontsize}, linewidths=.5,
                     cbar_kws=cbar_param)
    plt.title(title, fontsize=20)
    plt.xlabel('Localisation', fontsize=12)
    plt.ylabel('Superfamily', fontsize=12)

    cbar = hm.collections[0].colorbar
    cbar.ax.tick_params(labelsize=9)

    hm.tick_params(axis='y', labelsize=ylabel_fontsize)
    plt.show()


# heatmap(superfam_de_005[['upstream', 'downstream', 'intron', 'cds']], 'P4 MITEs (for DE genes < 0.05)', 10, 10, {"pad": 0.02, "shrink": 0.5})
heatmap(superfam_all_genes[['upstream', 'downstream', 'intron', 'cds']], 'P4 MITEs (for all genes)', 10, 10, {"pad": 0.02, "shrink": 0.5})
# heatmap(families_de_005[['upstream', 'downstream', 'intron', 'cds']], 'P4 MITEs (for DE genes < 0.05)',
#     7, 7, {"pad": 0.01, "shrink": 0.35})
# heatmap(families_all_genes[['upstream', 'downstream', 'intron', 'cds']], 'P4 MITEs (for all genes)',
#         5,5, {"pad": 0.01, "shrink": 0.35})
