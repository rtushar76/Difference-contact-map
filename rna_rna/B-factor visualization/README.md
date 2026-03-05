# Mapping Δ-Contact Values to RNA Structure (B-factor Visualization)

## Overview

The RNA--RNA contact difference matrix generated in this project
contains pairwise Δ-contact values for all nucleotide pairs. For large
RNAs such as the **16S rRNA**, this matrix becomes extremely large and
difficult to interpret directly as a heatmap.

This script converts the pairwise Δ-contact matrix into **per-nucleotide
cumulative Δ values** and writes them into the **B-factor column of a
PDB file**. This allows the contact differences to be visualized
directly on the RNA structure using molecular graphics software such as
**PyMOL, ChimeraX, or VMD**.

Each nucleotide receives a **single scalar value** representing the sum
of all Δ-contacts involving that nucleotide, which is then assigned to
the B-factor field of all atoms in that residue.

This enables intuitive visualization of regions of the RNA that **gain
or lose contacts between the two compared structures**.

------------------------------------------------------------------------

## Input Files

### Δ-Contact Matrix

    delta_contact_matrix.csv

Generated from the **RNA--RNA contact difference script**.

Matrix format example:

            A:1   A:2   A:3
    A:1      0     5    -2
    A:2      5     0     1
    A:3     -2     1     0

Rows and columns correspond to **nucleotide identifiers formatted as**:

    Chain:ResidueNumber

Example:

    A:150

------------------------------------------------------------------------

### RNA Structure

    8g7r_16srrna_chain_fixed.pdb

A reference **16S rRNA structure** used as the structural scaffold for
visualization.

The structure must use **chain IDs and residue numbering consistent with
the Δ-contact matrix**.

------------------------------------------------------------------------

## Output File

    rna_with_deltaB_vcmut-vcst.pdb

A modified PDB file where:

-   each nucleotide is assigned a **B-factor equal to its cumulative
    Δ-contact value**
-   all atoms within a nucleotide receive the **same B-factor value**

This file can be loaded into visualization software to display contact
differences on the RNA structure.

------------------------------------------------------------------------

## How the Script Works

### 1. Load the Δ-contact matrix

The script reads the CSV matrix containing pairwise contact differences.

### 2. Compute cumulative Δ per nucleotide

For each nucleotide:

    Δ_total(i) = Σ Δ(i,j) for all j ≠ i

This produces a **single scalar value per nucleotide** representing the
overall contact gain or loss.

### 3. Map values to the RNA structure

Each nucleotide in the input PDB is identified using:

    ChainID:ResidueNumber

Example:

    A:1453

The cumulative Δ value is assigned to the **B-factor field of all atoms
in that residue**.

### 4. Write the modified structure

The output PDB retains the original coordinates but contains **updated
B-factors representing Δ-contact values**.

------------------------------------------------------------------------

## Visualization

The resulting structure can be visualized in **PyMOL** or **ChimeraX**.

Example in PyMOL:

``` python
load rna_with_deltaB_vcmut-vcst.pdb
spectrum b, blue_white_red
```

Suggested interpretation:

  Color   Meaning
  ------- ------------------
  Blue    gain of contacts
  White   no change
  Red     loss of contacts

------------------------------------------------------------------------

## Why This Step Is Necessary

The **RNA--RNA Δ-contact matrix for the 16S rRNA is extremely large
(\~1500 × 1500)**. Direct heatmap visualization often obscures
structural patterns.

By projecting the matrix into **per-nucleotide Δ values**, this script
enables:

-   rapid identification of **regions with major contact
    rearrangements**
-   **structural visualization** of RNA remodeling
-   intuitive mapping of changes onto the ribosome architecture

------------------------------------------------------------------------

## Dependencies

Required Python packages:

    pandas
    biopython

Install with:

    pip install pandas biopython

------------------------------------------------------------------------

## Usage

Run the script directly:

    python map_delta_contacts_to_bfactor.py

The script will generate:

    rna_with_deltaB_vcmut-vcst.pdb

------------------------------------------------------------------------

## Notes

-   Residues not present in the Δ-contact matrix receive a **default
    B-factor of 0**.
-   The script assumes **consistent chain identifiers and residue
    numbering** between the matrix and the PDB file.
