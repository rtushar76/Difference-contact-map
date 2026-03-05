import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec

# Load downsampled matrix
df = pd.read_csv("delta_contact_matrix.csv", index_col=0)

# Clean up index and columns
df.index = df.index.astype(str).str.extract(r'(\d+)')[0].astype(int)
df.columns = df.columns.astype(str).str.extract(r'(\d+)')[0].astype(int)
df = df.sort_index().sort_index(axis=1)

# Prepare data for plotting
X, Y = np.meshgrid(df.columns, df.index)
Z = df.values

# Mutations to annotate
mutations = {426: "U426C (426)", 906: "A906G (906)", 1460: "A1460C (1459)", 1488: "G1488A (1487)", 1526: "C1526A (1525)"}

# Set up figure layout
fig = plt.figure(figsize=(14, 12))
gs = gridspec.GridSpec(1, 3, width_ratios=[0.96, 0.01, 0.03], wspace=0.40)

# Main plot axis
ax = fig.add_subplot(gs[0])

# Colorbar axis
cax = fig.add_subplot(gs[2])

# Keep the heatmap square
ax.set_aspect("equal")

# Define vmax
vmax = 30

# Plot heatmap using seismic_r so blue=positive and red=negative
pc = ax.pcolormesh(
    X,
    Y,
    Z,
    cmap=plt.cm.seismic_r,
    vmin=-vmax,
    vmax=vmax,
    shading='auto'
)

# Colorbar
cbar = fig.colorbar(pc, cax=cax)
cbar.set_label("ΔContacts", fontsize=18, labelpad=-17)
cbar.ax.tick_params(labelsize=16)

# Move the colorbar axis slightly to the right
pos = cax.get_position()
cax.set_position([pos.x0 + 0.02, pos.y0, pos.width, pos.height])

# Set sparse ticks (every 150 nucleotides)
x_ticks = np.arange(df.columns.min(), df.columns.max() + 1, 150)
y_ticks = np.arange(df.index.min(), df.index.max() + 1, 150)
ax.set_xticks(x_ticks)
ax.set_yticks(y_ticks)
ax.set_xticklabels([str(x) for x in x_ticks], rotation=90, fontsize=16)
ax.set_yticklabels([str(y) for y in y_ticks], fontsize=16)

# Axis labels
ax.set_xlabel("16S rRNA Nucleotide Index", fontsize=18)
ax.set_ylabel("16S rRNA Nucleotide Index", fontsize=18)
ax.set_title("ΔContacts within 16S rRNA (VC S4.4 - VC ST)", fontsize=20)

# Gridlines
ax.grid(which='major', color='black', linestyle=':', linewidth=0.4)

# Mutation labels and dashed lines
offsets = {
    426: +5,
    906: +5,
    1460: +5,
    1488: +5,
    1526: +5
}
y_shifts = {
    1460: -14,
    1526: +10
}

for nt, label in mutations.items():
    ax.axhline(nt, color='black', linestyle='--', linewidth=1)

    y_shift = y_shifts.get(nt, 0)
    x_offset = 10  # reasonable spacing on the right

    ax.text(
        df.columns.max() + x_offset,
        nt + y_shift,
        label,
        va='center',
        ha='left',
        fontsize=14
    )

# Final layout
plt.subplots_adjust(
    left=0.2,
    right=0.85,
    top=0.9,
    bottom=0.2
)

# Save figure
plt.savefig("intra_rna_delta_contacts_aggregated_ecmut_ecwt.png", dpi=600)
plt.show()