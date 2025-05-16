import numpy as np
import pandas as pd
from skbio.stats.ordination import pcoa
from scipy.spatial.distance import pdist, cdist

df = pd.read_csv('../files/heatmaps/cm_data/P1_all_genes_with_MITEs_cm_general_data_final_1.csv',sep='\t')
dmatrix = df.iloc[:, 1:-3].to_numpy()
print(type(dmatrix), dmatrix)

pdistance = cdist(dmatrix, dmatrix, 'euclidean')
print("P-distance", '\n', pdistance)

pcoa_result = pcoa(pdistance, number_of_dimensions=2)

