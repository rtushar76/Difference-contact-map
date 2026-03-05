# RNA–RNA (Intra-rRNA) Contact Difference Analysis

This script computes **intra-RNA (RNA–RNA) contact differences** between two 16S rRNA structures by comparing a **test structure (VC-MUT)** against a **reference structure (VC-ST)**.

The output is a **Δ contact matrix** (test − reference) that quantifies how **RNA–RNA contacts within the rRNA** change in the mutant relative to the reference.

---

## Overview

The workflow performs the following steps:

1. Parse a **Clustal alignment** to build nucleotide correspondence between VC-MUT and VC-ST
2. Load modeled residue indices from both RNA PDBs
3. Load per-nucleotide **Q-scores** for both structures
4. Build a nucleotide mapping (VC → reference) and filter mapped positions by:
   - nucleotide modeled in both structures
   - Q-score ≥ threshold in both structures
5. Compute **intra-RNA contact matrices** independently for each structure (atom–atom contacts summed per residue pair)
6. Symmetrize each matrix
7. Align both matrices using the VC numbering and compute:

**ΔC = C(VC-MUT) − C(VC-ST)**

8. Write outputs (delta matrix + residue mapping + exclusion log)

---

## Structures Compared

**Test structure (VC-MUT)**  
- RNA PDB: `vcmut_16srrna.pdb`  
- Q-scores: `vcs44_16srrna.txt`

**Reference structure (VC-ST)**  
- RNA PDB: `vcst_16srrna.pdb`  
- Q-scores: `vcst_16srrna.txt`

Interpretation of Δ values:

- **Positive** → contact gained/increased in mutant  
- **Negative** → contact lost/decreased in mutant  
- **Zero** → no change  

---

## Input Files

### 1) Clustal alignment

`clustal.aln`

Assumption used by the script:

- Alignment sequence 1 (index 0) = **reference** numbering (VC-ST / EC-style numbering in the script variable names)
- Alignment sequence 2 (index 1) = **test** numbering (VC-MUT)

Only aligned positions without gaps are mapped (VC index → reference index).

The script prints:

- number of aligned nucleotide pairs parsed from the alignment

---

### 2) RNA PDBs

- **Test (VC-MUT):** `vcmut_16srrna.pdb`
- **Reference (VC-ST):** `vcst_16srrna.pdb`

The script loads all atoms from all chains present in the PDBs, but the alignment-based mapping is constructed on nucleotide indices and Q-scores are queried as chain **A** by default (see Assumptions below).

---

### 3) Q-score files

- **Test (VC-MUT):** `vcs44_16srrna.txt`
- **Reference (VC-ST):** `vcst_16srrna.txt`

Expected line format (CryoSPARC-style):

```
Chain  ResNum  ...  Qavg
A      1405    ...  0.67
```

The script reads:
- `chain = parts[0]`
- `resnum = int(parts[1])`
- `q = float(parts[3])`

Residues are retained only if **both** structures satisfy:

- `Q ≥ 0.35` (default)

---

## Contact Definition (RNA–RNA)

The script computes an **intra-RNA contact matrix** by iterating over all RNA atoms.

Two atoms are considered in contact if their distance satisfies:

- `0.0 Å ≤ distance ≤ 4.0 Å`  (default)

Matrix meaning:

- **Rows / columns:** RNA residues labeled as `Chain:ResNum` (e.g., `A:1405`)
- **Values:** number of atom–atom contacts counted between the residue pair

Notes:
- Self-contacts (`A:1405` vs `A:1405`) are skipped.
- After computation, the matrix is **symmetrized** by `M = M + M.T`.

---

## Filtering Rules

A mapped nucleotide (VC idx ↔ reference idx) is excluded if:

- Either nucleotide is **unmodeled** in its respective RNA PDB (“Unmodeled”)
- Either nucleotide has **Q-score below threshold** (“Low Q-score”)

All excluded mapping pairs are logged.

---

## Output Files

### 1) Delta intra-RNA contact matrix

`delta_contact_matrix.csv`

This is a square matrix with:

- **Rows:** retained RNA residues (VC numbering, labeled `A:<VC_resnum>`)
- **Columns:** same as rows
- **Values:** Δ contact counts

Entry meaning:

`delta[i, j] = (atom-contact count between residues i and j in VC-MUT) − (atom-contact count between residues i and j in VC-ST)`

---

### 2) Residue mapping

`residue_mapping.txt`

Tab-separated mapping used to align matrices:

```
VC_residue    EC_residue
1405          1415
...
```

(Here “EC_residue” is the reference numbering used in the script variable names; in this context it corresponds to VC-ST.)

---

### 3) Excluded nucleotides log

`excluded_nucleotides.txt`

Tab-separated with columns:

- `VC_residue`
- `EC_residue` (reference residue index)
- `Reason` (Unmodeled / Low Q-score)

---

## Parameters (Defaults)

- `qscore_threshold = 0.35`
- `lower_cutoff = 0.0 Å`
- `upper_cutoff = 4.0 Å`

These are set in the script and can be adjusted there.

---

## Dependencies

Python packages required:

- `numpy`
- `pandas`
- `biopython`

Install via:

```
pip install numpy pandas biopython
```

---

## Running the Script

Place all required input files in the expected paths and run:

```
python rna_rna_difference_contact_map.py
```

The script prints:

- parsed alignment pair count
- number of nucleotides retained after Q filtering
- progress messages for matrix computation
- file saving messages

---

## Notes / Assumptions

1. VC-MUT is the **test** structure; VC-ST is the **reference** structure.
2. Alignment order is assumed to be:
   - sequence 0 = reference (VC-ST)
   - sequence 1 = test (VC-MUT)
3. Q-scores are queried using chain ID **A** (`("A", resnum)`).
4. The script checks whether nucleotides are modeled by searching residue indices across **all chains** present in each RNA PDB.
5. Output is an **intra-RNA atom-contact** delta matrix aligned into **VC numbering**.
