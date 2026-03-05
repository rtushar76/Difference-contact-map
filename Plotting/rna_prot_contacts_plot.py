import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
from matplotlib.ticker import FixedLocator, FixedFormatter

# Load CSV data
df = pd.read_csv("delta_contact_matrix.csv", index_col=0)

# Melt and process
df_long = df.stack().reset_index()
df_long.columns = ['RNA', 'Protein', 'Delta']
df_long = df_long[df_long['Delta'] != 0]
df_long['NucResi'] = df_long['RNA'].apply(lambda x: int(x.split(":")[1]))

# Map protein chains to names
chain_to_protein = {
    "B": "uS2", "C": "uS3", "D": "uS4", "E": "uS5", "F": "bS6", "G": "uS7", "H": "uS8",
    "I": "uS9", "J": "uS10", "K": "uS11", "L": "uS12", "M": "uS13", "N": "uS14",
    "O": "uS15", "P": "bS16", "Q": "uS17", "R": "bS18", "S": "uS19", "T": "bS20", "U": "uS21"
}
df_long['ProteinName'] = df_long['Protein'].apply(lambda x: chain_to_protein.get(x.split(":")[0], "unknown"))

# Aggregate delta contacts
agg = df_long.groupby(['NucResi', 'ProteinName'])['Delta'].sum().unstack(fill_value=0)

# Axis labels
xlabels = ["uS2", "uS3", "uS4", "uS5", "bS6", "uS7", "uS8", "uS9", "uS10", "uS11", "uS12", "uS13",
           "uS14", "uS15", "bS16", "uS17", "bS18", "uS19", "bS20", "uS21"]
ylabels = sorted(agg.index)
agg = agg.reindex(columns=xlabels, fill_value=0)
data = agg.loc[ylabels, xlabels].values

# Create figure and gridspec identical to intra-RNA
fig = plt.figure(figsize=(14, 12))
gs = gridspec.GridSpec(1, 3, width_ratios=[0.96, 0.01, 0.03], wspace=0.3)

# Main heatmap axis
ax = fig.add_subplot(gs[0])

# Colorbar axis on the right
cax = fig.add_subplot(gs[2])

# Keep the heatmap square
#ax.set_aspect("equal")

# Heatmap
cmap = plt.cm.seismic_r
vmax = np.max(np.abs(data))
X, Y = np.meshgrid(np.arange(len(xlabels)+1), np.array(ylabels + [ylabels[-1]+1]))
pc = ax.pcolormesh(
    X - 0.5,
    Y - 0.5,
    data,
    cmap=cmap,
    vmin=-vmax,
    vmax=vmax,
    shading='flat'
)

# Colorbar
cbar = fig.colorbar(pc, cax=cax)
cbar.set_label("ΔContacts", fontsize=18, labelpad=-17)
cbar.ax.tick_params(labelsize=16)

# X ticks and labels
ax.xaxis.set_major_locator(FixedLocator(np.arange(len(xlabels))))
ax.xaxis.set_major_formatter(FixedFormatter(xlabels))
plt.setp(ax.get_xticklabels(), rotation=90, fontsize=16)

# Y ticks and labels
yticks = list(range(min(ylabels), max(ylabels) + 1, 150))
ax.yaxis.set_major_locator(FixedLocator(yticks))
ax.yaxis.set_major_formatter(FixedFormatter([str(y) for y in yticks]))
plt.setp(ax.get_yticklabels(), fontsize=16)

# Labels and title
ax.set_xlabel("Ribosomal Proteins", fontsize=18)
ax.set_ylabel("16S rRNA Nucleotide Index", fontsize=18)
ax.set_title("ΔContacts per 16S rRNA Nucleotide vs r-Protein (VC S4.4 - VC ST)", fontsize=20, pad=25)

# Grid
ax.grid(True, linestyle=':', linewidth=0.3, zorder=1)

# Mutation labels and dashed lines
mutations = {
    426: "U426C (426)",
    906: "A906G (906)",
    1460: "A1460C (1459)",
    1488: "G1488A (1487)",
    1526: "C1526A (1525)"
}

offsets = {
    426: -0.25,
    906: -0.25,
    1460: -0.25,
    1488: -0.25,
    1526: -0.25
}
y_shifts = {
    1460: -10,
    1526: +8
}

for nt, label in mutations.items():
    ax.axhline(nt, color='black', linestyle='--', linewidth=1)

    y_shift = y_shifts.get(nt, 0)
    x_offset = offsets.get(nt, +5)

    ax.text(
        len(xlabels) + x_offset,
        nt + y_shift,
        label,
        va='center',
        ha='left',
        fontsize=14
    )

# Final layout
plt.subplots_adjust(
    left=0.2,
    right=0.93,
    top=0.9,
    bottom=0.2
)

# Save figure
plt.savefig("rna_protein_delta_contacts.png", dpi=600)
plt.show()
