import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.colors as mcolors
import matplotlib.ticker as ticker
from matplotlib.patches import Rectangle


# 1. Input files
mites_bed = "../files/mite_density_data/EL10_mites_to_calc_density.bed"           # format: chr start end
chrom_lengths_file = "../files/mite_density_data/chrom_lengths.csv"  # format: chr length
bin_size = 100000                 # 100 kb windows

# 2. Load data
mites = pd.read_csv(mites_bed, sep="\t", header=None, names=["chr", "start", "end"])
chrom_lengths = pd.read_csv(chrom_lengths_file, sep="\t", header=None, names=["chr", "length"])

# ensure numeric
mites["start"] = mites["start"].astype(int)
mites["end"] = mites["end"].astype(int)
chrom_lengths["length"] = chrom_lengths["length"].astype(int)

# 3. Compute MITE density per bin
density_data = []

for _, row in chrom_lengths.iterrows():
    chrom = row["chr"]
    chrom_len = row["length"]
    n_bins = int(np.ceil(chrom_len / bin_size))

    # select mites on this chromosome
    subset = mites[mites["chr"] == chrom]

    # create bins
    bins = np.arange(0, chrom_len + bin_size, bin_size)
    mite_counts, _ = np.histogram(subset["start"], bins=bins)

    # compute density as number per Mb (optional)
    density = mite_counts / (bin_size / 100000)  # per 100kb

    density_data.append(pd.DataFrame({
        "chr": chrom,
        "bin_start": bins[:-1],
        "bin_end": bins[1:],
        "density": density
    }))

density_df = pd.concat(density_data, ignore_index=True)

# 4. Plot as chromosome stripes
chrom_order = chrom_lengths["chr"].tolist()
density_df["chr"] = pd.Categorical(density_df["chr"], categories=chrom_order, ordered=True)
density_df = density_df.sort_values(["chr", "bin_start"])

fig_height = max(6, len(chrom_lengths) * 0.8)  # auto-adjust plot height
fig, ax = plt.subplots(figsize=(14, fig_height))

# Normalize color scale
norm = mcolors.Normalize(vmin=0, vmax=density_df["density"].max())

# first version of stripes - without outline
# stripe_height = 0.6
#
# for i, chrom in enumerate(chrom_order):
#     chrom_data = density_df[density_df["chr"] == chrom]
#     chrom_len = chrom_lengths.loc[chrom_lengths["chr"] == chrom, "length"].values[0]
#     n_bins = len(chrom_data)
#     ax.barh(
#         y=[i]*len(chrom_data),
#         width=chrom_data["bin_end"] - chrom_data["bin_start"],
#         left=chrom_data["bin_start"],
#         color=plt.cm.BuPu(norm(chrom_data["density"])),
#         height=0.8,
#         edgecolor='none'
#     )
# # Outline around each chromosome
#     ax.add_patch(Rectangle(
#         (0, i - stripe_height/2),
#         chrom_len,
#         stripe_height,
#         fill=False,
#         color="black",
#         linewidth=0.8
#     ))

# nice version of chromosomes - with outline adjusted properly
stripe_height = 0.8  # match the bar height

for i, chrom in enumerate(chrom_order):
    chrom_data = density_df[density_df["chr"] == chrom]
    chrom_len = chrom_lengths.loc[chrom_lengths["chr"] == chrom, "length"].values[0]

    # Draw chromosome density stripes
    ax.barh(
        y=[i]*len(chrom_data),
        width=chrom_data["bin_end"] - chrom_data["bin_start"],
        left=chrom_data["bin_start"],
        color=plt.cm.BuPu(norm(chrom_data["density"])),
        height=stripe_height,
        edgecolor='none'
    )

    # Draw outline *around* the whole chromosome
    ax.add_patch(Rectangle(
        (0, i - stripe_height/2),  # shift slightly lower
        chrom_len,
        stripe_height,  # make it a bit taller
        fill=False,
        edgecolor="black",
        linewidth=1.0,
        zorder=3  # ensure it draws on top
    ))

# Set x-axis ticks every 10 Mb
max_length = chrom_lengths["length"].max()
xticks = np.arange(0, max_length + 1, 10_000_000)   # every 10 Mb
xticklabels = [f"{x/1_000_000:.0f}" for x in xticks]
ax.set_xticklabels(xticklabels, fontsize=18)

# Disable scientific notation on x-axis
ax.xaxis.set_major_formatter(ticker.FuncFormatter(lambda x, pos: f"{int(x/1_000_000)}"))

ax.set_yticks(range(len(chrom_order)))
ax.set_yticklabels(chrom_order, fontsize=18)
ax.set_xlabel("Genomic position [Mb]", fontsize=25)
ax.set_ylabel("Chromosome", fontsize=25)
ax.invert_yaxis()  # optional, top to bottom order
plt.title("EL10 - MITE density per chromosome", fontsize=30)

# Add colorbar
sm = plt.cm.ScalarMappable(cmap='BuPu', norm=norm)
sm.set_array([])
cbar = plt.colorbar(sm, ax=ax, orientation="vertical", pad=0.02)
cbar.set_label("MITE density (per 100 kb)", fontsize=18)
ax.set_xlim(0, max_length)
plt.tight_layout()
plt.show()
