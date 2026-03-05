RNA–Protein Contact Difference Analysis

This repository contains a Python script to compute RNA–protein contact
differences between two ribosome structures.

The script compares a test structure (VC-MUT) against a reference
structure (VC-ST) and quantifies contact gains and losses between rRNA
nucleotides and ribosomal protein residues.

------------------------------------------------------------------------

Overview

The workflow performs the following steps:

1.  Parse a Clustal alignment to build nucleotide correspondence between
    the two structures
2.  Filter residues using Q-scores and modeled residue checks
3.  Compute RNA–protein contact matrices for both structures
4.  Align matrices using canonical RNA numbering
5.  Compute the Δ contact matrix (test − reference)

The resulting matrix describes how RNA–protein contacts change in the
mutant relative to the reference.

------------------------------------------------------------------------

Structures Compared

The analysis compares two structures:

Test structure (VC-MUT) Mutant ribosome structure being analyzed.

Reference structure (VC-ST) Baseline ribosome structure used for
comparison.

The contact difference is defined as:

ΔC = C(VC-MUT) − C(VC-ST)

Interpretation:

Positive → Contact gained in mutant
Negative → Contact lost in mutant
Zero → No change

------------------------------------------------------------------------

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

------------------------------------------------------------------------

2. PDB Files

Test Structure (VC-MUT)

vcmut_16srrna.pdb
vcmut_rproteins.pdb

Reference Structure (VC-ST)

vcst_16srrna.pdb
vcst_rproteins.pdb

------------------------------------------------------------------------

3. Q-Score Files

VC-MUT

vcs44_16srrna.txt
vcs44_rproteins.txt

VC-ST

vcst_16srrna.txt
vcst_rproteins.txt

Expected format:

Chain Residue … Qscore

A 1405 … 0.67
A 1406 … 0.55

Residues are retained only if:

Q-score ≥ 0.35

------------------------------------------------------------------------

Contact Definition

Two atoms are considered in contact if the distance satisfies:

0 Å ≤ distance ≤ 4 Å

All atom-atom contacts between each RNA nucleotide and protein residue
are counted.

Matrix meaning:

Rows → RNA residues
Columns → Protein residues
Values → Number of atom–atom contacts

------------------------------------------------------------------------

Filtering Steps

RNA residues are excluded if:

• Not modeled in either structure
• Q-score below threshold in either structure
• Alignment gap

Logged in:

excluded_nucleotides.txt

Protein residues are excluded if:

Q-score < threshold in either structure

Logged in:

excluded_proteins.txt

------------------------------------------------------------------------

Matrix Alignment

After computing contact matrices for both structures:

VC-MUT matrix
VC-ST matrix

They are aligned using the RNA mapping derived from the Clustal
alignment.

Protein residues are restricted to the intersection of residues present
in both matrices.

------------------------------------------------------------------------

Output Files

delta_contact_matrix.csv

Rows → RNA residues (VC numbering)
Columns → Protein residues
Values → Δ contact counts

Interpretation:

Δ > 0 → Contact gained in mutant
Δ < 0 → Contact lost in mutant
Δ = 0 → No change

------------------------------------------------------------------------

residue_mapping.txt

VC_residue EC_residue 1405 1415 1406 1416

------------------------------------------------------------------------

excluded_nucleotides.txt

Example:

VC_residue EC_residue Reason 1403 1413 Low Q-score 1407 1417 Unmodeled

------------------------------------------------------------------------

excluded_proteins.txt

Example:

Chain ResNum Reason D 53 Low Q-score E 108 Low Q-score

------------------------------------------------------------------------

Dependencies

numpy
pandas
biopython

Install with:

pip install numpy pandas biopython

------------------------------------------------------------------------

Running the Script

Place all input files in the same directory and run:

python contact_difference_analysis.py

The script prints progress messages indicating:

• parsed alignment positions
• retained residues after filtering
• excluded residues
• contact matrix calculation

------------------------------------------------------------------------

Expected Use Cases

This analysis can be used to:

• identify RNA–protein interaction changes
• study mutational adaptation of ribosomes
• quantify contact gains and losses across rRNA

The resulting Δ matrix can be used for:

• B-factor coloring in PyMOL
• contact heatmaps
• per-residue aggregation analyses

------------------------------------------------------------------------

Assumptions

The script assumes:

1.  RNA chain ID is A
2.  Clustal alignment order is EC first, VC second
3.  Residue numbering between PDB and Q-score files is consistent
4.  Q-score files follow a CryoSPARC-style format
