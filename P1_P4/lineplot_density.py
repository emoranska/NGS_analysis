import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt

sns.set_theme(style="ticks", palette="deep")
mites_dens = pd.read_csv('../files/mite_density_data/P4_mite_density_data.csv', sep='\t')

# dens = sns.lineplot(data=mites_dens, x="range", y="all_mites_dens")
mites_dens_wide = mites_dens.pivot(index="range", columns="chr", values="all_mites_dens")
print(mites_dens_wide.to_string(max_rows=10))

plt.figure(figsize=(20, 10))

dens_wide = sns.lineplot(data=mites_dens_wide)
plt.title("P4", fontsize=35)
plt.xlabel("Chromosome length [Mb]", fontsize=23)
plt.ylabel("MITE density [per 100 kb]", fontsize=23)
dens_wide.tick_params(axis='y', labelsize=20)
dens_wide.tick_params(axis='x', labelsize=20)
plt.legend(fontsize=14, loc='lower right')
plt.show()
