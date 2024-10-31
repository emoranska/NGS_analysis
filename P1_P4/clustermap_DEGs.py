import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns

data_stow_3 = pd.read_csv('../files/Stow3_clustermap_data.csv', sep='\t', header=0, index_col=0)


def clustermap_degs(data, mite_fam):
    genes = data.index
    sample_names = data.columns

    d_st3 = data[data.columns[0:6]]

    # to check how to add color labels to rows according to values from 'mite_loc'
    mite_loc_labels = data['mite_loc']
    mite_loc_pal = sns.color_palette("Spectral_r", 3)
    mite_loc_lut = dict(zip(mite_loc_labels.unique(), mite_loc_pal))
    row_colors = mite_loc_labels.map(mite_loc_lut)

    degs_clustermap = sns.clustermap(d_st3, cmap="bwr", xticklabels=sample_names, yticklabels=genes, z_score=0, center=0,
                                     col_cluster=False, figsize=(8,9), row_colors=row_colors, cbar_kws={'label': 'row Z-score'},
                                     cbar_pos=(.99, .65, .02, .15))

    ax = degs_clustermap.ax_heatmap
    ax.set_ylabel("")

    ax.set_title(f"{mite_fam} MITEs associated with DEGs", fontsize=20)

    for label in mite_loc_labels.unique():
        degs_clustermap.ax_col_dendrogram.bar(0, 0, color=mite_loc_lut[label],
                                label=label, linewidth=0)
    degs_clustermap.ax_col_dendrogram.legend(loc="center left", ncol=4)


# check for changing y label colours
    genes_labels = data['ins_de']
    genes_lut = dict(zip(genes_labels.unique(), 'bwr'))

    for tick_label in degs_clustermap.ax_heatmap.axes.get_yticklabels():
        tick_text = tick_label.get_text()
        gene_name = genes_labels.loc[tick_text]
        tick_label.set_color(genes_lut[gene_name])

    # degs_clustermap.savefig(f'../files/clustermaps/{mite_fam}.pdf')
    plt.show()

clustermap_degs(data_stow_3, "Stow3")
