import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import matplotlib.gridspec as gridspec

# Map chain IDs to ribosomal protein names
chain_to_protein = {
    "B": "uS2", "C": "uS3", "D": "uS4", "E": "uS5", "F": "bS6", "G": "uS7", "H": "uS8",
    "I": "uS9", "J": "uS10", "K": "uS11", "L": "uS12", "M": "uS13", "N": "uS14",
    "O": "uS15", "P": "bS16", "Q": "uS17", "R": "bS18", "S": "uS19", "T": "bS20", "U": "uS21"
}

# Load the matrix CSV
df = pd.read_csv("delta_contact_matrix.csv", index_col=0)

# Rename the rows and columns to protein names
df.rename(index=chain_to_protein, columns=chain_to_protein, inplace=True)

# Define consistent order
ordered_labels = [
    "uS2", "uS3", "uS4", "uS5", "bS6", "uS7", "uS8", "uS9", "uS10", "uS11",
    "uS12", "uS13", "uS14", "uS15", "bS16", "uS17", "bS18", "uS19", "bS20", "uS21"
]
df = df.reindex(index=ordered_labels, columns=ordered_labels, fill_value=0)

# Extract data
data = df.values
xlabels = df.columns.tolist()
ylabels = df.index.tolist()

# Figure and gridspec (same layout as RNA script)
fig = plt.figure(figsize=(14, 12))
gs = gridspec.GridSpec(1, 3, width_ratios=[0.96, 0.01, 0.03], wspace=0.3)

# Colorbar axis
ax = fig.add_subplot(gs[0])

# Main plot axis
cax = fig.add_subplot(gs[2])

# Keep the heatmap square
ax.set_aspect("equal")

# Color scale
cmap = plt.cm.seismic_r
vmax = np.max(np.abs(data))
if vmax == 0:
    vmax = 1

X, Y = np.meshgrid(np.arange(len(xlabels)+1), np.arange(len(ylabels)+1))

# Heatmap
pc = ax.pcolormesh(
    X - 0.5, Y - 0.5, data,
    cmap=cmap,
    vmin=-vmax,
    vmax=vmax,
    shading="flat"
)

# Colorbar
cbar = fig.colorbar(pc, cax=cax)
cbar.set_label("ΔContacts", fontsize=18, labelpad=-25)
cbar.ax.tick_params(labelsize=16)

# Ticks
from matplotlib.ticker import FixedLocator, FixedFormatter

ax.xaxis.set_major_locator(FixedLocator(np.arange(len(xlabels))))
ax.xaxis.set_major_formatter(FixedFormatter(xlabels))
ax.yaxis.set_major_locator(FixedLocator(np.arange(len(ylabels))))
ax.yaxis.set_major_formatter(FixedFormatter(ylabels))

plt.setp(ax.get_xticklabels(), rotation=90, fontsize=16)
plt.setp(ax.get_yticklabels(), fontsize=16)

# Labels and title
ax.set_xlabel("Ribosomal Proteins", fontsize=18)
ax.set_ylabel("Ribosomal Proteins", fontsize=18)
ax.set_title("ΔContacts between Ribosomal Proteins (VC S4.4 - VC ST)", fontsize=20)

# Grid
ax.grid(True, linestyle=":", linewidth=0.3, zorder=1)

# Layout and save
plt.subplots_adjust(
    left=0.3,
    right=0.85,
    top=0.9,
    bottom=0.2
)

plt.savefig("intra_protein_delta_contacts_q_filtered.png", dpi=600)
plt.show()