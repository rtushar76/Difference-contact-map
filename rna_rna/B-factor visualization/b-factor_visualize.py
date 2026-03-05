import pandas as pd
from Bio.PDB import PDBParser, PDBIO
from collections import defaultdict

# Inputs
delta_csv = "delta_contact_matrix.csv"
input_pdb = "8g7r_16srrna_chain_fixed.pdb"
output_pdb = "rna_with_deltaB_vcmut-vcst.pdb"

# Load Δ-contact matrix
df = pd.read_csv(delta_csv, index_col=0)

# Initialize delta map
delta_map = defaultdict(float)

# Parse the matrix and sum delta values for each nucleotide (excluding diagonal)
for res1 in df.index:
    for res2 in df.columns:
        if res1 != res2:
            delta = df.loc[res1, res2]
            delta_map[res1] += delta
            delta_map[res2] += delta

# Load PDB structure
parser = PDBParser(QUIET=True)
structure = parser.get_structure("rna", input_pdb)

# Apply B-factors per residue
for model in structure:
    for chain in model:
        for residue in chain:
            resid_str = f"{chain.id}:{residue.id[1]}"
            delta_value = delta_map.get(resid_str, 0.0)
            for atom in residue:
                atom.bfactor = delta_value

# Save output
io = PDBIO()
io.set_structure(structure)
io.save(output_pdb)
print(f"Saved: {output_pdb}")
