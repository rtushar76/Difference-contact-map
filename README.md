## License

This project is licensed under the MIT License.  
See the LICENSE file for details.

# Ribosome Contact Difference Analysis

This repository contains scripts used to compute and visualize structural contact differences between ribosomal structures.

The workflow compares a **test structure (VC-MUT)** with a **reference structure (VC-ST)** and quantifies changes in molecular contacts across three interaction classes:

• RNA–RNA contacts within 16S rRNA  
• RNA–protein contacts between rRNA and ribosomal proteins  
• Protein–protein contacts between ribosomal proteins  

Each analysis produces a Δ contact matrix:

ΔC = C_test − C_reference

Interpretation:

Positive → Contact gained in mutant  
Negative → Contact lost in mutant  
Zero → No change

The resulting matrices can be visualized using the plotting scripts included in the repository.

------------------------------------------------------------
Repository Structure
------------------------------------------------------------

rna_rna/
    Scripts for computing contact differences within 16S rRNA.

rna_protein/
    Scripts for computing RNA–protein contact differences.

protein_protein/
    Scripts for computing ribosomal protein contact differences.

plotting/
    Scripts used to visualize Δ contact matrices as heatmaps.

Each directory contains its own README explaining the specific workflow.

------------------------------------------------------------
Dependencies
------------------------------------------------------------

Required Python packages:

numpy  
pandas  
biopython  
matplotlib  
scipy  

Install with:

pip install numpy pandas biopython matplotlib scipy

------------------------------------------------------------
License
------------------------------------------------------------

This project is licensed under the MIT License.
See the LICENSE file for details.

------------------------------------------------------------
Author
------------------------------------------------------------

Tushar Raskar  
Structural Biology / Cryo-EM
