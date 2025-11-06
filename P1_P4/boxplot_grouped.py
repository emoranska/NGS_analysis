import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt

sns.set_theme(style="ticks", palette="deep")
f, ax = plt.subplots(figsize=(14, 10))
ax.set_yscale("log")

# mites_plants = pd.read_csv('../files/boxplots/P1_MITEs_in_genes_plants_boxplot_data.csv', sep='\t')
mites_plants = pd.read_csv('../files/boxplots/P1_P4_uc_plants_to_boxplot.csv', sep='\t')

# plt.figure(figsize=(14, 10))
# sns.boxplot(x="mite_group", y="Mite_sum", hue="mite_loc", data=mites_plants)
# sns.despine(offset=10, trim=True)
bplot = sns.boxplot(x="mite_loc", y="mite_sum", data=mites_plants, color="white", fliersize=0, whis=(0, 100), width=.6,
                    boxprops=dict(linewidth=2, edgecolor="black"), medianprops=dict(linewidth=1.5, color="black"),
                    whiskerprops=dict(linewidth=1.5, color="black"), capprops=dict(linewidth=1.5, color="black"),
                    native_scale=True)

# Overlay scatter (or swarm) points colored by F2_family
strplot = sns.stripplot(
    x="mite_loc",
    y="mite_sum",
    data=mites_plants,
    hue="F2_family",
    dodge=False,        # ensures all points stay within one box
    jitter=0.25,        # spread the points
    palette="deep",     # choose a nice color palette
    #alpha=0.8,
    size=9,
    linewidth=1,
    edgecolor="black"
    )

# bplot = sns.boxplot(x="mite_loc", y="mite_sum", hue="F2_family", data=mites_plants)
# sns.violinplot(x="mite_loc", y="mite_sum", hue="pop", data=mites_plants, split=True, inner="quart")

plt.title("uc", fontsize=35)
plt.xlabel(None)
plt.ylabel("Number of MITEs (log scale)", fontsize=28)
strplot.tick_params(axis='y', labelsize=25)
strplot.tick_params(axis='x', labelsize=28, rotation=15)
plt.legend(fontsize=25, loc="lower left")
plt.show()
