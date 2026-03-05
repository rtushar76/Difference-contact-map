RNA–Protein Contact Difference Analysis

This repository contains a Python script to compute RNA–protein contact differences between two ribosome structures.

The script compares a test structure (VC-MUT) against a reference structure (VC-ST) and quantifies contact gains and losses between rRNA nucleotides and ribosomal protein residues.

Overview

The workflow performs the following steps:

Parse a Clustal alignment to build nucleotide correspondence between the two structures

Filter residues using Q-scores and modeled residue checks

Compute RNA–protein contact matrices for both structures

Align matrices using canonical RNA numbering

Compute the Δ contact matrix (test − reference)

The resulting matrix describes how RNA–protein contacts change in the mutant relative to the reference.

Structures Compared

The analysis compares two structures:

Test structure (VC-MUT)
Mutant ribosome structure being analyzed.

Reference structure (VC-ST)
Baseline ribosome structure used for comparison.

The contact difference is defined as:

ΔC = C(VC-MUT) − C(VC-ST)

Interpretation:

Value	Meaning
Positive	Contact gained in mutant
Negative	Contact lost in mutant
Zero	No change
Input Files
1. Sequence Alignment
clustal.aln

Clustal alignment containing the 16S rRNA sequences of both structures.

Assumptions used by the script:

Sequence 1 → EC / reference numbering

Sequence 2 → VC / mutant numbering

Only positions aligned without gaps are used to build residue mapping.

Example mapping:

VC nucleotide 1405 ↔ EC nucleotide 1415
2. PDB Files
Test Structure (VC-MUT)
vcmut_16srrna.pdb
vcmut_rproteins.pdb

Contains:

RNA atoms

ribosomal protein atoms

Reference Structure (VC-ST)
vcst_16srrna.pdb
vcst_rproteins.pdb

Contains:

RNA atoms

ribosomal protein atoms

3. Q-Score Files

Residue-level Q-scores for both structures.

VC-MUT
vcs44_16srrna.txt
vcs44_rproteins.txt
VC-ST
vcst_16srrna.txt
vcst_rproteins.txt

Expected format:

Chain  Residue  ...  Qscore
A      1405     ...  0.67
A      1406     ...  0.55

Residues are retained only if:

Q-score ≥ 0.35

This threshold is defined in the script and can be modified.

Contact Definition

Two atoms are considered in contact if the distance satisfies:

0 Å ≤ distance ≤ 4 Å

All atom-atom contacts between each RNA nucleotide and protein residue are counted.

The resulting matrix has:

Dimension	Meaning
Rows	RNA residues
Columns	Protein residues
Values	Number of atom–atom contacts
Filtering Steps
RNA residues

A nucleotide is excluded if:

Not modeled in either structure

Q-score below threshold in either structure

Alignment gap

Excluded residues are recorded in:

excluded_nucleotides.txt
Protein residues

A protein residue is excluded if:

Q-score < threshold in either structure

Excluded residues are logged in:

excluded_proteins.txt
Matrix Alignment

After computing contact matrices for both structures:

VC-MUT matrix
VC-ST matrix

They are aligned using the RNA mapping derived from the Clustal alignment.

Protein residues are restricted to the intersection of residues present in both matrices.

Output Files
1. Delta Contact Matrix
delta_contact_matrix.csv

Matrix dimensions:

Rows	RNA residues (VC numbering)
Columns	Protein residues
Values	Δ contact counts

Interpretation:

Value	Meaning
Δ > 0	Contact gained in mutant
Δ < 0	Contact lost in mutant
Δ = 0	No change
2. Residue Mapping
residue_mapping.txt

Mapping between canonical residues:

VC_residue    EC_residue
1405          1415
1406          1416
...
3. Excluded Nucleotides
excluded_nucleotides.txt

Lists residues removed due to:

alignment gaps

missing structural model

low Q-score

Example:

VC_residue EC_residue Reason
1403       1413       Low Q-score
1407       1417       Unmodeled
4. Excluded Protein Residues
excluded_proteins.txt

Example:

Chain ResNum Reason
D     53     Low Q-score
E     108    Low Q-score
Dependencies

Python packages required:

numpy
pandas
biopython

Install via:

pip install numpy pandas biopython
Running the Script

Place all input files in the same directory and run:

python contact_difference_analysis.py

The script will print progress messages indicating:

parsed alignment positions

retained residues after filtering

excluded residues

contact matrix calculation

Expected Use Cases

This analysis can be used to:

identify RNA–protein interaction changes

study mutational adaptation of ribosomes

quantify contact gains and losses across rRNA

The resulting Δ matrix can be used for:

B-factor coloring in PyMOL

contact heatmaps

per-residue aggregation analyses

Assumptions

The script assumes:

RNA chain ID is A

Clustal alignment order is:

EC first
VC second

Residue numbering between PDB and Q-score files is consistent

Q-score files follow a CryoSPARC-style format
