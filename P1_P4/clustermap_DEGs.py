import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns

data_stow_3 = pd.read_csv('../files/heatmaps/cm_data/Stow3_cm_data.csv', sep='\t',
                          header=0, index_col=0)
data_stow_36 = pd.read_csv('../files/heatmaps/cm_data/Stow36_cm_data.csv', sep='\t', header=0,
                           index_col=0)
data_stow_40 = pd.read_csv('../files/heatmaps/cm_data/Stow40_cm_data.csv', sep='\t', header=0,
                           index_col=0)
data_tourist_1 = pd.read_csv('../files/heatmaps/cm_data/Tourist1_cm_data.csv', sep='\t', header=0,
                             index_col=0)
data_tourist_17 = pd.read_csv('../files/heatmaps/cm_data/Tourist17_cm_data.csv', sep='\t', header=0,
                              index_col=0)
data_uc_33 = pd.read_csv('../files/heatmaps/cm_data/uc33_cm_data.csv', sep='\t', header=0,
                         index_col=0)


def clustermap_degs(data, mite_fam):
    genes = data.index
    sample_names = data.columns

    data_tpm = data[data.columns[0:6]]

    # variables for adding color labels to rows according to values from 'mite_loc'
    mite_loc_labels = data['mite_loc']
    mite_loc_options = ['upstream', '5UTR', 'intron', 'cds', '3UTR', 'downstream']
    mite_loc_pal = sns.color_palette("Spectral_r",7)
    # mite_loc_lut = dict(zip(mite_loc_labels.unique(), mite_loc_pal))
    mite_loc_lut = dict(zip(mite_loc_options, mite_loc_pal))
    print(mite_loc_lut)
    row_colors = mite_loc_labels.map(mite_loc_lut)

    # create clustermap
    degs_clustermap = sns.clustermap(data_tpm, cmap="bwr", xticklabels=sample_names, yticklabels=genes, z_score=0,
                                     center=0, col_cluster=False, figsize=(6, 7), row_colors=row_colors,
                                     cbar_kws={'label': 'row Z-score'}, cbar_pos=(.95, .65, .02, .15))

    ax = degs_clustermap.ax_heatmap

    # remove label from y-axis
    ax.set_ylabel("")

    # set main title and axis labels font size
    ax.set_title(f"{mite_fam} MITEs associated with DEGs", fontsize=18)
    degs_clustermap.tick_params(axis='x', labelsize=8)
    degs_clustermap.tick_params(axis='y', labelsize=8)

    # add color labels to rows according to values from 'mite_loc'
    for label in mite_loc_labels.unique():
        degs_clustermap.ax_col_dendrogram.bar(0, 0, color=mite_loc_lut[label], label=label, linewidth=0)
    degs_clustermap.ax_col_dendrogram.legend(loc="center", ncol=4)

    # changing y-axis label colours according to 'ins_de' values
    genes_labels = data['ins_de']
    print(type(genes_labels), genes_labels)
    genes_lut = dict(zip(genes_labels.unique(), 'rbg'))
    print(type(genes_lut), print(genes_lut))

    for tick_label in degs_clustermap.ax_heatmap.axes.get_yticklabels():
        tick_text = tick_label.get_text()
        print(type(tick_text), tick_text)
        gene_name = genes_labels.loc[tick_text]
        print(type(gene_name), gene_name)
        tick_label.set_color(genes_lut[gene_name])

    degs_clustermap.savefig(f'../files/heatmaps/clustermaps/{mite_fam}_no_repeats.pdf')
    plt.show()


# clustermap_degs(data_stow_3, "Stow3")
# clustermap_degs(data_stow_36, "Stow36")
# clustermap_degs(data_stow_40, "Stow40")
# clustermap_degs(data_tourist_1, "Tourist1")
# clustermap_degs(data_tourist_17, "Tourist17")
clustermap_degs(data_uc_33, "uc33")
