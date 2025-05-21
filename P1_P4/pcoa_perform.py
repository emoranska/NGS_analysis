import pandas as pd
from skbio.stats.ordination import pcoa
from scipy.spatial.distance import pdist, cdist
import matplotlib.pyplot as plt
import seaborn as sns

# df = pd.read_csv('../files/heatmaps/cm_data/P1_all_genes_with_MITEs_cm_general_data_final_1.csv',sep='\t')
# df = pd.read_csv('../files/heatmaps/cm_data/P4_all_genes_with_MITEs_cm_general_data_final_1.csv',sep='\t')
df = pd.read_csv('../files/heatmaps/cm_data/P4_all_genes_cm_general_data_final.csv', sep='\t')
dmatrix = df.iloc[:, 1:-3].to_numpy()
print(type(dmatrix), dmatrix)

pdistance = cdist(dmatrix, dmatrix, 'euclidean')
print("P-distance", '\n', pdistance)

pcoa_result = pcoa(pdistance, number_of_dimensions=2)
coordinates = pcoa_result.samples
print('Coordinates:', '\n', coordinates)

df_pcoa = coordinates[['PC1', 'PC2']]
df_pcoa['mite_loc'] = df['mite_loc'].to_numpy()
df_pcoa = df_pcoa.set_index('mite_loc')
df_pcoa['ins_de'] = df['ins_de'].to_numpy()
print('\n\n', df_pcoa)

# fig, ax = plt.subplots()
# df_pcoa.plot('PC1', 'PC2', kind='scatter', ax=ax) #c='ins_de', colormap='viridis')
# plt.title('PCoA Plot')
#
# for k, v in df_pcoa.iterrows():
#     ax.annotate(k, v)

# color_dict = dict({'no_DE': 'green',
#                   'up': 'red',
#                   'down': 'blue'})

g = sns.scatterplot(x="PC1", y="PC2", hue='mite_loc', style='ins_de',
              data=df_pcoa, #palette=color_dict,
                   legend='full')

plt.title('P4 - all genes')
plt.show()
