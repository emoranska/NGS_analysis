import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns

superfam_de_005 = (pd.read_csv('../files/P4_DE_005_superfam_counts.csv', sep='\t',
                               dtype={'MITE_sum': np.int32, 'upstream': np.int32, 'downstream': np.int32,
                                      'intron': np.int32, 'cds': np.int32}).set_index('superfamily'))

superfam_de_005_data = superfam_de_005[['upstream', 'downstream', 'intron', 'cds']]
print(superfam_de_005_data.to_string())

superfam_all_genes = (pd.read_csv('../files/P4_all_genes_superfam_counts.csv', sep='\t').
                      set_index('superfamily'))
superfam_all_genes_data = superfam_all_genes[['upstream', 'downstream', 'intron', 'cds']]


def heatmap(data, title):
    # działa
    # df.style.background_gradient(cmap='Blues').set_properties(**{'font-size': '20px'}).to_excel('styled.xlsx')

    # plt.pcolor(df)
    # # plt.yticks(np.arange(0.5, len(rpk_data.index), 1), rpk_data.index)
    # # plt.xticks(np.arange(0.5, len(rpk_data.columns), 1), rpk_data.columns)
    # plt.show()

    # # Displaying dataframe as a heatmap
    # # with diverging colourmap as RdYlBu
    # plt.imshow(df, cmap="RdYlBu")
    #
    # # Displaying a color bar to understand
    # # which color represents which range of data
    # plt.colorbar()
    #
    # # # Assigning labels of x-axis
    # # # according to dataframe
    # # plt.xticks(range(len(df)), df.columns)
    # #
    # # # Assigning labels of y-axis
    # # # according to dataframe
    # # plt.yticks(range(len(df)), df.index)
    #
    # # Displaying the figure
    # plt.show()

    fig, ax = plt.subplots(figsize=(12, 7))
    sns.heatmap(data, cmap="YlGnBu", annot=True, fmt="0.0f", square=True)
    plt.title(title, fontsize=20)
    plt.xlabel('Localisation', fontsize=12)
    plt.ylabel('Superfamily', fontsize=12)
    plt.show()


# heatmap(superfam_de_005_data, 'P4 MITEs (for DE genes < 0.05)')
heatmap(superfam_all_genes_data, 'P4 MITEs (for all genes)')
