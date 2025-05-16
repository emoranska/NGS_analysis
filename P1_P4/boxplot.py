import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns

superfam_all_mites = pd.read_csv('../files/boxplots/P1_boxplot_data.csv', sep='\t')
superfam_de005 = pd.read_csv('../files/boxplots/P4_DE005_boxplot_data.csv', sep='\t')


def boxplot(data):
    plt.figure(figsize=(10, 10))
    bplot = sns.boxplot(data, x="MITE_group", y="MITE_sum")
    bplot.set_yticks(range(len(data)-50))

    bplot.tick_params(axis='y', labelsize=18)
    plt.ylabel("Number of copies per MITE family", fontsize=20)
    bplot.tick_params(axis='x', labelsize=18)
    plt.xlabel(None)
    plt.ylim(0, 20)
    plt.show()


# boxplot(superfam_all_mites)
boxplot(superfam_de005)
