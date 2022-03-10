# PDownB

The RCSB PDB was searched for crystal structures to sample. Crystal structures were used instead of NMR models to simplify the data analysis process; NMR structures often contain multiple models that could be problematic to process.

The criteria for PDB sampling were:

- Resolution > 2 Å
- RWork between 0-0.20
- RFree between 0-0.25
- Chain length 30 - 300 residues
- No ligands
- Minimum 2 alpha helices
- Remove redundancy by using representative structures at 70 % sequence identity (i.e. use one structure to represent all lysozymes).

The search produced **2367 entries** for which the coordinate (pdb) files were downloaded for analysis. 
- In "PDB" you will find these PDB files as they were submitted to the RCSB PDB.
- In "CA" you will find the alpha-carbon backbone for each protein, as a .txt file. If the protein has multiple chains, we selected the first available chain.
- In "XYZ" you will find the coordinate files for these backbone curves of each protein.
- In "SS" you will find the secondary structure classification for each protein, where "-" denotes linker, "H" denotes alpha-helix, and "S" denotes beta-sheets.
- In "Writhes" you will find the writhe data of each protein. The format of this is subsection length ---------- writhe. We start from length 5 (the smallest possible to calculate the writhe) and increase in steps of 3.
