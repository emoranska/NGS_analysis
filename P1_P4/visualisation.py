import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns

df_superfam = (pd.read_csv('../files/superfam_counts.csv', sep='\t',
                 dtype={'MITE_sum': np.int32, 'upstream': np.int32, 'downstream': np.int32, 'intron': np.int32,
                        'cds': np.int32, 'de_results_up': np.int32, 'de_results_down': np.int32}).
      set_index('superfamily'))

df = df_superfam[['upstream', 'downstream', 'intron', 'cds']]

print(df.to_string())

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
sns.heatmap(df, cmap="YlGnBu", annot=True, fmt="0.0f", square=True)
plt.title('P4 MITEs (for DE genes < 0.05)', fontsize=20)
plt.xlabel('Localisation', fontsize=12)
plt.ylabel('Superfamily', fontsize=12)
plt.show()
