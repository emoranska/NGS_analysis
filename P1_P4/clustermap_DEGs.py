import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns

data_stow_3 = pd.read_csv('../files/Stow3_clustermap_data.csv', sep='\t', header=0, index_col=0)

genes = data_stow_3.index
sample_names = data_stow_3.columns

degs_clustermap = sns.clustermap(data_stow_3, xticklabels=sample_names, yticklabels=genes, z_score=0, col_cluster=False,
                                 figsize=(7,9))

plt.show()
