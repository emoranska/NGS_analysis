import pandas as pd
import numpy as np
from tqdm import tqdm
import matplotlib.pyplot as plt

# -----------------------------
# 1. Load input data
# -----------------------------
mites = pd.read_csv("../files/permutation/P1_mites_in_bins_to_permutation_no_dupl.bed", sep="\t",
                    header=None, names=["chr", "start", "end"])
genes = pd.read_csv("../files/permutation/P1_genes_in_bins_to_permutation.bed", sep="\t",
                    header=None, names=["chr", "start", "end", "gene_id", "is_DEG"])

mites['start'] = mites['start'].astype(int)
mites['end'] = mites['end'].astype(int)
genes['start'] = genes['start'].astype(int)
genes['end'] = genes['end'].astype(int)
# -----------------------------
# 2. Define association criterion
# -----------------------------
# Here: a gene is "associated" if a MITE overlaps within ±2 kb of its region
WINDOW = 2000


def count_mite_associated_genes(mites, genes):
    """Count how many DEGs are associated with MITEs."""
    assoc = []
    for _, gene in genes.iterrows():
        subset = mites[mites["chr"] == gene["chr"]]
        overlap = subset[
            (subset["end"] >= gene["start"] - WINDOW) &
            (subset["start"] <= gene["end"] + WINDOW)
        ]
        if not overlap.empty:
            assoc.append(gene["is_DEG"])
    if len(assoc) == 0:
        return 0
    return sum(assoc)  # number of DEGs associated with MITEs


# Observed number of MITE–DEG associations
observed = count_mite_associated_genes(mites, genes)
print(f"Observed MITE–DEG associations: {observed}")

# -----------------------------
# 3. Permutation test
# -----------------------------
n_permutations = 1000

np.random.seed(42)

# Store null distribution
null_counts = []

for i in tqdm(range(n_permutations), desc="Permuting"):
    shuffled_genes = genes.copy()
    shuffled_genes["is_DEG"] = np.random.permutation(shuffled_genes["is_DEG"].values)
    count = count_mite_associated_genes(mites, shuffled_genes)
    null_counts.append(count)

null_counts = np.array(null_counts)
np.savetxt("../files/permutation/P1_null_counts.csv", null_counts, delimiter="\t")

# -----------------------------
# 4. Compute empirical p-value
# -----------------------------
p_value = (np.sum(null_counts >= observed) + 1) / (n_permutations + 1)
print(f"Empirical p-value: {p_value:.4g}")

# -----------------------------
# 5. Optional: visualize
# -----------------------------

# correct number of MITEs associated wit DEGs according to data in Supplementary Materials 2
observed = 36
p_value = (np.sum(null_counts >= observed) + 1) / (n_permutations + 1)
print(f"Empirical p-value with correct number of MITEs associated with DEGs: {p_value}", type(p_value))
# print(f"Empirical p-value: {p_value:.4g}")

plt.hist(null_counts, bins=30, color="gray", alpha=0.7)
plt.axvline(observed, color="red", linestyle="dashed", linewidth=2, label="Observed")
plt.xlabel("Number of MITE–DEG associations (randomized)")
plt.ylabel("Frequency")
plt.title("Permutation Test for MITE–DEG Enrichment")
plt.legend(loc='upper center')
plt.show()
