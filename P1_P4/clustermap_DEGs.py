import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
from scipy.cluster import hierarchy

# data_stow_3 = pd.read_csv('../files/heatmaps/cm_data/Stow3_cm_data_correct.csv', sep='\t',
#                           header=0, index_col=0)
# data_stow_36 = pd.read_csv('../files/heatmaps/cm_data/Stow36_cm_data_correct.csv', sep='\t', header=0,
#                            index_col=0)
# data_stow_40 = pd.read_csv('../files/heatmaps/cm_data/Stow40_cm_data_correct.csv', sep='\t', header=0,
#                            index_col=0)
# data_tourist_1 = pd.read_csv('../files/heatmaps/cm_data/Tourist1_cm_data_correct.csv', sep='\t', header=0,
#                              index_col=0)
# data_tourist_17 = pd.read_csv('../files/heatmaps/cm_data/Tourist17_cm_data_correct.csv', sep='\t', header=0,
#                               index_col=0)
# data_uc_33 = pd.read_csv('../files/heatmaps/cm_data/uc33_cm_data_correct.csv', sep='\t', header=0,
#                          index_col=0)
# data_deg_with_mites = pd.read_csv('../files/heatmaps/cm_data/P1_DE005_with_MITEs_cm_general_data_final_1.csv',
#                                   sep='\t', header=0, index_col=0)
data_all_genes_with_mites = pd.read_csv('../files/heatmaps/cm_data/'
                                        'P1_all_genes_with_MITEs_cm_general_data_final_1.csv',
                                        sep='\t', header=0, index_col=0)
# data_all_degs = pd.read_csv('../files/heatmaps/cm_data/P1_all_DEGs_cm_general_data_final.csv',
#                             sep='\t', header=0, index_col=0)
# data_all_genes = pd.read_csv('../files/heatmaps/cm_data/P4_all_genes_cm_general_data_final.csv',
#                              sep='\t', header=0, index_col=0)


def clustermap_degs(data, mite_fam):
    genes = data.index
    sample_names = data.columns

    data_tpm = data[data.columns[0:6]]

    # variables for adding color labels to rows according to values from 'mite_loc'
    mite_loc_labels = data['mite_loc']
    # mite_loc_options = ['upstream', '5UTR', 'intron', 'cds', '3UTR', 'downstream']
    mite_loc_options = ['upstream', '5UTR', 'intron', 'cds', '3UTR', 'downstream', 'no_TE']
    mite_loc_pal = sns.color_palette("Spectral_r", 7)
    mite_loc_lut = dict(zip(mite_loc_options, mite_loc_pal))
    print(mite_loc_lut)
    row_colors = mite_loc_labels.map(mite_loc_lut)

    # create clustermap
    # degs_clustermap = sns.clustermap(data_tpm, cmap="bwr", xticklabels=sample_names, yticklabels=genes, z_score=0,
    #                                  center=0, col_cluster=False, figsize=(6, 7), row_colors=row_colors,
    #                                  cbar_kws={'label': 'row Z-score'}, cbar_pos=(.82, .65, .02, .15))


    degs_clustermap = sns.clustermap(data_tpm, cmap="bwr", xticklabels=sample_names, yticklabels=genes, z_score=0,
                                     center=0, col_cluster=False, figsize=(6, 7), row_colors=row_colors,
                                     cbar_kws={'label': 'row Z-score'}, cbar_pos=(.82, .65, .02, .15),
                                     method='complete', metric='cityblock')

    ax = degs_clustermap.ax_heatmap

    # remove label from y-axis
    ax.set_ylabel("")

    # set main title and axis labels font size
    # ax.set_title(f"{mite_fam} MITEs associated with DEGs", fontsize=18)
    ax.set_title(f"{mite_fam} MITEs associated with all genes", fontsize=18)
    # ax.set_title(f"{mite_fam} DEGs", fontsize=18)
    # ax.set_title(f"{mite_fam} genes", fontsize=18)
    degs_clustermap.tick_params(axis='x', labelsize=8)
    degs_clustermap.tick_params(axis='y', labelsize=2)

    # add color labels to rows according to values from 'mite_loc'
    for label in mite_loc_labels.unique():
        degs_clustermap.ax_col_dendrogram.bar(0, 0, color=mite_loc_lut[label], label=label, linewidth=0)
    degs_clustermap.ax_col_dendrogram.legend(loc="center", ncol=4)

    # changing y-axis label colours according to 'ins_de' or 'de_all' values
    genes_labels = data['ins_de']
    # genes_labels = data['de_all']
    print(type(genes_labels), genes_labels)
    # genes_labels_options = ['ins_up', 'ins_down']
    # genes_labels_options = ['up', 'down']
    genes_labels_options = ['up', 'down', 'no_DE']
    # genes_labels_options = ['up', 'down', 'TE_no_DE', 'no_TE']
    # genes_lut = dict(zip(genes_labels.unique(), 'rbg'))
    # genes_lut = dict(zip(genes_labels_options, 'rbgy'))
    genes_lut = dict(zip(genes_labels_options, 'rbg'))
    print(type(genes_lut), print(genes_lut))

    for tick_label in degs_clustermap.ax_heatmap.axes.get_yticklabels():
        tick_text = tick_label.get_text()
        # print(type(tick_text), tick_text)
        gene_name = genes_labels.loc[tick_text]
        print(type(gene_name), gene_name)
        tick_label.set_color(genes_lut[gene_name])

    # degs_clustermap.savefig(f'../files/heatmaps/clustermaps/{mite_fam}_genes.pdf')
    plt.show()


# clustermap_degs(data_stow_3, "Stow3")
# clustermap_degs(data_stow_36, "Stow36")
# clustermap_degs(data_stow_40, "Stow40")
# clustermap_degs(data_tourist_1, "Tourist1")
# clustermap_degs(data_tourist_17, "Tourist17")
# clustermap_degs(data_uc_33, "uc33")
# clustermap_degs(data_deg_with_mites, ' P1_all')
clustermap_degs(data_all_genes_with_mites, ' P1_all')
# clustermap_degs(data_all_degs, 'P1_all')
# clustermap_degs(data_all_genes, 'P4_all')
