import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns

superfam_all_mites = pd.read_csv('../files/boxplots/P4_boxplot_mean_supfam_data.csv', sep='\t')
superfam_de005 = pd.read_csv('../files/boxplots/P4_DE005_boxplot_data.csv', sep='\t')
supfam_plants_all_mites = pd.read_csv('../files/boxplots/P1_boxplot_plants_all_mites_supfam_data.csv',
                                      sep='\t')
supfam_all_vs_plants = pd.read_csv('../files/boxplots/P4_boxplot_all_mites_vs_plants_supfam_data_1.csv',
                                   sep='\t')


def boxplot(data):
    plt.figure(figsize=(10, 10))
    bplot = sns.boxplot(data, x="MITE_group", y="MITE_sum")
    # bplot = sns.boxplot(data, x="MITE_group", y="MITE_all_vs_MITE_plant")
    # bplot.set_yticks(range(len(data)-50))

    bplot.tick_params(axis='y', labelsize=18)
    # plt.ylabel("Number of copies per MITE family", fontsize=20)
    plt.ylabel("Number of copies per genome", fontsize=23)
    # plt.ylabel("MITE_all / MITE_plant", fontsize=23)
    bplot.tick_params(axis='x', labelsize=23)
    plt.xlabel(None)
    # plt.ylim(0, 20)
    # plt.ylim(0, 2500)
    # plt.ylim(1, 1.5)
    plt.show()


# boxplot(superfam_all_mites)
# boxplot(superfam_de005)
boxplot(supfam_plants_all_mites)
# boxplot(supfam_all_vs_plants)
