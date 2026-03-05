# --- USER PARAMETERS ---

# Input files
clustal_alignment_file = "clustal.aln"

vc_rna_pdb = "vcmut_16srrna.pdb"
ec_rna_pdb = "vcst_16srrna.pdb"

qscore_vc_rna_file = "vcs44_16srrna.txt"
qscore_ec_rna_file = "vcst_16srrna.txt"

# Q-score cutoff
qscore_threshold = 0.35

# Contact parameters
lower_cutoff = 0.0
upper_cutoff = 4.0

# --- FUNCTIONS ---

def parse_alignment(alignment_file):
    alignment = AlignIO.read(alignment_file, "clustal")
    seq_vc = alignment[1].seq
    seq_ec = alignment[0].seq
    mapping = {}
    vc_idx = 0
    ec_idx = 0
    for i in range(len(seq_vc)):
        vc_char = seq_vc[i]
        ec_char = seq_ec[i]
        if vc_char != "-":
            vc_idx += 1
        if ec_char != "-":
            ec_idx += 1
        if vc_char != "-" and ec_char != "-":
            mapping[vc_idx] = ec_idx
    print(f"Parsed {len(mapping)} aligned nucleotide pairs from Clustal alignment.")
    return mapping

def get_pdb_residues(pdb_file):
    parser = PDBParser(QUIET=True)
    structure = parser.get_structure("pdb", pdb_file)
    residues = defaultdict(set)
    for model in structure:
        for chain in model:
            for res in chain:
                residues[chain.id].add(res.id[1])
    return residues

def load_qscores(path):
    qdict = {}
    with open(path) as f:
        for line in f:
            if line.startswith("#") or not line.strip():
                continue
            parts = line.strip().split()
            chain = parts[0]
            resnum = int(parts[1])
            q = float(parts[3])
            qdict[(chain, resnum)] = q
    return qdict

def load_atoms(pdb_file):
    parser = PDBParser(QUIET=True)
    structure = parser.get_structure("pdb", pdb_file)
    atoms = []
    for model in structure:
        for chain in model:
            for res in chain:
                for atom in res:
                    atoms.append({
                        'chain': chain.id,
                        'resi': res.id[1],
                        'coord': atom.coord
                    })
    return atoms

def compute_intra_rna_contact_matrix(rna_atoms, excluded_rna):
    residues = sorted(set((a['chain'], a['resi']) for a in rna_atoms if (a['chain'], a['resi']) not in excluded_rna))
    matrix = pd.DataFrame(0, index=[f"{c}:{r}" for c,r in residues],
                             columns=[f"{c}:{r}" for c,r in residues])
    for atom1 in rna_atoms:
        res1 = (atom1['chain'], atom1['resi'])
        if res1 in excluded_rna:
            continue
        label1 = f"{atom1['chain']}:{atom1['resi']}"
        for atom2 in rna_atoms:
            res2 = (atom2['chain'], atom2['resi'])
            if res2 in excluded_rna:
                continue
            label2 = f"{atom2['chain']}:{atom2['resi']}"
            if label1 == label2:
                continue  # Skip self-contacts
            dist = np.linalg.norm(atom1['coord'] - atom2['coord'])
            if lower_cutoff <= dist <= upper_cutoff:
                matrix.at[label1, label2] += 1
    return matrix

# --- MAIN EXECUTION ---

print("Parsing canonical alignment...")
alignment_mapping = parse_alignment(clustal_alignment_file)

print("Loading residues modeled in PDBs...")
vc_rna_residues = get_pdb_residues(vc_rna_pdb)
ec_rna_residues = get_pdb_residues(ec_rna_pdb)

print("Loading Q-scores...")
q_vc_rna = load_qscores(qscore_vc_rna_file)
q_ec_rna = load_qscores(qscore_ec_rna_file)

# Map canonical numbering to PDB numbering, with logging of exclusions
rna_mapping = {}
rna_exclusions = []

for vc_idx, ec_idx in alignment_mapping.items():
    modeled_vc = any(vc_idx in rset for rset in vc_rna_residues.values())
    modeled_ec = any(ec_idx in rset for rset in ec_rna_residues.values())

    if not modeled_vc or not modeled_ec:
        reason = "Unmodeled"
        rna_exclusions.append((vc_idx, ec_idx, reason))
        continue

    q_vc = q_vc_rna.get(("A", vc_idx), 0)
    q_ec = q_ec_rna.get(("A", ec_idx), 0)

    if q_vc < qscore_threshold or q_ec < qscore_threshold:
        reason = "Low Q-score"
        rna_exclusions.append((vc_idx, ec_idx, reason))
        continue

    rna_mapping[vc_idx] = ec_idx

print(f"Total retained nucleotides after Q filtering: {len(rna_mapping)}")

print("Loading RNA atoms...")
rna_atoms_vc = load_atoms(vc_rna_pdb)
rna_atoms_ec = load_atoms(ec_rna_pdb)

print("Computing intra-RNA contact matrices...")
matrix_vc = compute_intra_rna_contact_matrix(rna_atoms_vc, excluded_rna=set())
matrix_ec = compute_intra_rna_contact_matrix(rna_atoms_ec, excluded_rna=set())

# Make matrices symmetric
matrix_vc = matrix_vc + matrix_vc.T
matrix_ec = matrix_ec + matrix_ec.T

# Align matrices by VC numbering, restrict to common nucleotides only
vc_labels = [f"A:{v}" for v in rna_mapping.keys()]
ec_labels = [f"A:{rna_mapping[v]}" for v in rna_mapping.keys()]

aligned_vc = matrix_vc.loc[vc_labels, vc_labels]
aligned_ec = matrix_ec.loc[ec_labels, ec_labels]
aligned_ec.index = vc_labels
aligned_ec.columns = vc_labels

# Compute delta
delta = aligned_vc - aligned_ec

print("Saving delta contact matrix...")
delta.to_csv("delta_contact_matrix.csv")

print("Saving residue mapping...")
with open("residue_mapping.txt", "w") as f:
    f.write("VC_residue\tEC_residue\n")
    for vc_idx in sorted(rna_mapping.keys()):
        f.write(f"{vc_idx}\t{rna_mapping[vc_idx]}\n")

print("Saving nucleotide exclusion log...")
with open("excluded_nucleotides.txt", "w") as f:
    f.write("VC_residue\tEC_residue\tReason\n")
    for v, e, reason in rna_exclusions:
        f.write(f"{v}\t{e}\t{reason}\n")

print("Done.")
