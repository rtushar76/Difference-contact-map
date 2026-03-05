# Protein–Protein (Chain–Chain) Contact Difference Analysis

This script computes **protein–protein contact differences** between two ribosome structures by comparing a **test structure (VC-MUT)** against a **reference structure (VC-ST)**.

Unlike an atom-contact map, this script produces a **chain–chain difference matrix** where each entry reflects the **difference in the number of residue–residue contacts** between two chains in the test vs reference structures.

---

## Overview

The workflow performs the following steps:

1. Load per-residue **Q-scores** for both structures
2. Load per-residue **atomic coordinates** from both protein PDBs
3. For each ribosomal protein chain, parse a **Clustal alignment** to build a residue index mapping (test ↔ reference)
4. Filter residue pairs by:
   - residue modeled in both structures
   - Q-score ≥ threshold in both structures
5. For each chain pair (chain1, chain2), compute a **residue–residue contact count** for each structure:
   - a residue pair counts as a contact if **any atom–atom distance < cutoff**
6. Compute the **Δ chain–chain contact matrix (test − reference)** and write outputs

---

## Structures Compared

**Test structure (VC-MUT)**  
- `vcmut_rproteins.pdb`  
- Q-scores: `vcs44_rproteins.txt`

**Reference structure (VC-ST)**  
- `vcst_rproteins.pdb`  
- Q-scores: `vcst_rproteins.txt`

The contact difference is defined as:

**ΔC = C(VC-MUT) − C(VC-ST)**

Interpretation:

- **Positive** → more residue–residue contacts in the mutant (gain)  
- **Negative** → fewer residue–residue contacts in the mutant (loss)  
- **Zero** → no net change  

---

## Input Files

### 1) Protein PDBs

- **Test (VC-MUT):** `vcmut_rproteins.pdb`  
- **Reference (VC-ST):** `vcst_rproteins.pdb`

Both PDBs must contain the **same set of chain IDs** (the script checks this and aborts if mismatched).

---

### 2) Q-score files

- **Test (VC-MUT):** `vcs44_rproteins.txt`  
- **Reference (VC-ST):** `vcst_rproteins.txt`

Expected per-line format (CryoSPARC-style):

```
Chain  ResNum  ...  Qavg
B      15      ...  0.62
```

The script reads:
- `chain = parts[0]`
- `resnum = int(parts[1])`
- `qavg = float(parts[3])`

Residues are retained only if **both** structures satisfy:

- `Qavg ≥ 0.4`  (default)

---

### 3) Per-chain Clustal alignments

A **separate Clustal alignment** is provided for each chain, e.g.

- `alignments/B.aln`
- `alignments/C.aln`
- ...
- `alignments/U.aln`

Assumption used by the script for each alignment file:

- Alignment sequence 1 = **test** chain numbering (VC-MUT)
- Alignment sequence 2 = **reference** chain numbering (VC-ST)

Only aligned positions without gaps are mapped.

---

## Contact Definition (Residue–Residue)

For a residue pair (r1 in chain1, r2 in chain2):

- Compute all atom–atom distances between the two residues using `scipy.spatial.distance.cdist`
- The pair is considered a **contact** if:

**min(atom–atom distance) < 4.0 Å** (default)

Important: the script counts **1 contact per residue pair** (not multiple contacts per atom pair).

---

## Filtering Rules

A residue mapping pair is excluded if:

- Either residue is **not modeled** in the corresponding PDB chain (“Unmodeled”)
- Either residue has **Qavg below threshold** in its respective structure (“Low Q-score”)

All exclusions are logged.

---

## Output Files

### 1) Delta chain–chain contact matrix

`delta_contact_matrix.csv`

This is a square matrix with:

- **Rows:** protein chain IDs (e.g., B–U)  
- **Columns:** protein chain IDs (e.g., B–U)  
- **Values:** Δ contact counts between the chain pair

Entry meaning:

- `matrix[ch1, ch2] = (# residue-pair contacts between ch1 and ch2 in VC-MUT) − (# residue-pair contacts between ch1 and ch2 in VC-ST)`

---

### 2) Excluded residue mapping log

`excluded_residues.txt`

Tab-separated with columns:

- `Chain`
- `Res1` (test residue index)
- `Res2` (reference residue index)
- `Reason` (Unmodeled / Low Q-score)

---

## Parameters (Defaults)

- `qscore_threshold = 0.4`
- `distance_cutoff = 4.0 Å`

These are set in the script and can be adjusted there.

---

## Dependencies

Python packages required:

- `numpy`
- `pandas`
- `biopython`
- `scipy`

Install via:

```
pip install numpy pandas biopython scipy
```

---

## Running the Script

Place all inputs in the expected paths and run:

```
python protein_protein_difference_map.py
```

The script prints:

- detected chains
- alignment parsing status
- contact matrix computation progress
- file saving messages

---

## Notes / Assumptions

1. VC-MUT is the **test** structure; VC-ST is the **reference** structure.
2. Chain IDs must match between the two PDBs.
3. Residue numbering in the Q-score files must match the PDB residue numbering.
4. Each chain alignment file is Clustal formatted and ordered as:
   - sequence 1 = test (VC-MUT)
   - sequence 2 = reference (VC-ST)
5. Output is **chain–chain** delta counts (not a residue-by-residue delta map).
