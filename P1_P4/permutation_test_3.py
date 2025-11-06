import pandas as pd
import numpy as np
from tqdm import tqdm
import matplotlib.pyplot as plt
import os

# 1. Load input data
mites = pd.read_csv("../files/permutation/P4_mites_in_bins_to_permutation_no_dupl.bed", sep="\t",
                    header=None, names=["chr", "start", "end"])
genes = pd.read_csv("../files/permutation/P4_genes_in_bins_to_permutation.bed", sep="\t",
                    header=None, names=["chr", "start", "end", "gene_id", "is_DEG"])

for df in (mites, genes):
    df["start"] = df["start"].astype(int)
    df["end"] = df["end"].astype(int)


# 2. Define parameters
WINDOW = 2000       # ±2 kb flanking window
N_PERMUTATIONS = 1000
OUTPUT_DIR = "../files/permutation/results/"
os.makedirs(OUTPUT_DIR, exist_ok=True)

np.random.seed(42)


# 3. Define association function
# Count DEGs and total genes associated with MITEs.
# A gene is 'associated' if any MITE overlaps or lies within ±window bp.
def count_mite_associated_genes(mites, genes, window=2000):
    deg_count = 0
    total_count = 0

    for _, gene in genes.iterrows():
        subset = mites[mites["chr"] == gene["chr"]]
        nearby = subset[
            (subset["end"] >= gene["start"] - window) &
            (subset["start"] <= gene["end"] + window)
        ]
        if not nearby.empty:
            total_count += 1
            if gene["is_DEG"] == 1:
                deg_count += 1

    return deg_count, total_count


# 4. Compute observed statistics
deg_near, total_near = count_mite_associated_genes(mites, genes, WINDOW)
observed_prop = deg_near / total_near if total_near > 0 else np.nan
print(f"Observed: {deg_near} DEGs / {total_near} genes near MITEs "
      f"(within ±{WINDOW} bp or overlapping) = {observed_prop:.4f}")


# 5. Permutation test
null_props = []

for i in tqdm(range(N_PERMUTATIONS), desc="Permuting"):
    shuffled = genes.copy()
    shuffled["is_DEG"] = np.random.permutation(shuffled["is_DEG"].values)
    deg_c, tot_c = count_mite_associated_genes(mites, shuffled, WINDOW)
    null_props.append(deg_c / tot_c if tot_c > 0 else np.nan)

null_props = np.array(null_props)
null_props = null_props[~np.isnan(null_props)]  # remove NaNs if any


# 6. Compute empirical p-value
p_value = (np.sum(null_props >= observed_prop) + 1) / (len(null_props) + 1)
print(f"Empirical p-value (combined overlap + ±{WINDOW} bp): {p_value:.4g}")

# 7. Save results
null_path = os.path.join(OUTPUT_DIR, f"P4_null_distribution_combined_{WINDOW}bp.csv")
summary_path = os.path.join(OUTPUT_DIR, f"P4_summary_combined_{WINDOW}bp.txt")

np.savetxt(null_path, null_props, delimiter="\t", fmt="%.6f")

with open(summary_path, "w") as f:
    f.write("MITE–DEG permutation test (combined overlap + flanking window)\n")
    f.write(f"Window: ±{WINDOW} bp\n")
    f.write(f"Permutations: {N_PERMUTATIONS}\n")
    f.write(f"Observed DEGs near MITEs: {deg_near}\n")
    f.write(f"Total genes near MITEs: {total_near}\n")
    f.write(f"Observed proportion: {observed_prop:.6f}\n")
    f.write(f"Empirical p-value: {p_value:.6f}\n")
    f.write(f"Null distribution file: {null_path}\n")

print(f"\nResults saved to:\n  - {null_path}\n  - {summary_path}")


# 8. Visualize results
null_props = pd.read_csv("../files/permutation/results/P4_null_distribution_combined_2000bp.csv")
observed_prop = 0.098

plt.figure(figsize=(12, 8))
plt.rc('xtick', labelsize=20)
plt.rc('ytick', labelsize=20)
plt.hist(null_props, bins=30, color="gray", alpha=0.7)
plt.axvline(observed_prop, color="red", linestyle="dashed", linewidth=2, label="Observed proportion")
plt.axvline(observed_prop, color="red", linestyle="dashed", linewidth=2, label="Observed proportion")
plt.xlabel("Proportion of DEGs associated with MITEs (randomized)", fontsize=25)
plt.ylabel("Frequency", fontsize=25)
# plt.title(f"Permutation Test: DEG enrichment near/overlapping MITEs (±{WINDOW} bp)")
plt.legend(loc='upper center')
plt.tight_layout()

plt.show()
