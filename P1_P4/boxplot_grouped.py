import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt

sns.set_theme(style="ticks", palette="deep")
# mites_plants = pd.read_csv('../files/boxplots/P1_MITEs_in_genes_plants_boxplot_data.csv', sep='\t')
mites_plants = pd.read_csv('../files/boxplots/P1_P4_all_MITEs_plants_to_boxplot.csv', sep='\t')

# sns.boxplot(x="mite_group", y="Mite_sum", hue="mite_loc", data=mites_plants)
# sns.despine(offset=10, trim=True)

plt.figure(figsize=(14, 10))
bplot = sns.boxplot(x="mite_loc", y="mite_sum", hue="F2_family", data=mites_plants)
# sns.violinplot(x="mite_loc", y="mite_sum", hue="pop", data=mites_plants, split=True, inner="quart")
plt.title("All MITEs", fontsize=35)
plt.xlabel(None)
plt.ylabel("Number of MITEs", fontsize=28)
bplot.tick_params(axis='y', labelsize=25)
bplot.tick_params(axis='x', labelsize=28, rotation=15)
plt.legend(fontsize=25)
plt.show()
