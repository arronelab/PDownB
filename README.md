# PDownB

The RCSB PDB was searched for crystal structures to sample. Crystal structures were used instead of NMR models to simplify the data analysis process; NMR structures often contain multiple models that could be problematic to process.

The criteria for PDB sampling were:

- Resolution > 2 Å
- Chain length 30 - 300 residues
- Remove redundancy by using representative structures at 70 % sequence identity (i.e. use one structure to represent all lysozymes).

The search produced **2367 entries** for which the coordinate (pdb) files were downloaded for analysis.
- **PDB**: the PDB files as they were submitted to the RCSB PDB for each protein.
- **CA**: the alpha-carbon backbone for each protein. If the protein has multiple chains, we selected the first available chain.
- **XYZ**: the xyz-coordinate files for the backbone curves of each protein.
- **SSE**: the secondary structure classification for each protein, where "-" denotes linker, "H" denotes alpha-helix, and "S" denotes beta-sheets.
- **Down**: the xyz-coordinate files for the downsampled proteins. The downsampling reduces the protein to just the endpoints of each secondary structure element.
- **DownSSE**: the secondary structure classification for each downsampled protein, where "-" denotes linker, "H" denotes alpha-helix, and "S" denotes beta-sheets.
- **Writhes**: the writhe data of each downsampled protein in the format: *Subsection Length* | *Writhe*. We start from length 5 (the smallest possible to calculate the writhe) and increase in steps of 3.

*Note: There are only 2333 entries in the downsampled case, as some proteins when downsampled were then too short to calculate a writhe for.*

- **curvature_torsion_calculator.py**: the script used to calculate the curvature-torsion distribution of the proteins in **PDB**. *Requires [pandas](https://pandas.pydata.org/) and [Biopython](https://biopython.org/).*
- **curvature_torsion_analysis.csv**: the output of **curvature_torsion_calculator.py** in the format: *Chain* | *Subsection* | *Curvature* | *Torsion*

It is possible to get distributions for each secondary structure element using [DSSP](https://swift.cmbi.umcn.nl/gv/dssp/). If this is something you are interested in then please get in touch for more details.

**c++Molecule** contains the code to make many predictions from one input file (currently human_SMARCAL1).
- You will first need to compile this code by running ```sh makeFileFinalPrediction.sh```. This should produce an executable called **predictStructure**.
- To run this executable, first make the config file by running ```python makePredictionFile.py```. You will be asked for three inputs:
  1. First, a project name. *e.g.* 2Changes4Times
  2. Next, the number of changes from the initial structure. This is the number of times the algorithm will run through the molecule and make changes to the whole structure. *e.g.* 2 as per the project name. *Note this is changes to the whole structure, that is if you choose to make just one change, it is each atom that will be changed once, not just one atom.* 
  3. Finally, the number of times you'd like to perform this. *e.g.* 4 as per the project name. *Note this is the number of times you go from the initial structure, and make the number of changes defined above. This is therefore the number of output files you will get.*
- You will now have a file **projectname_config.sh**. You are now ready to make predictions with ```sh projectname_config.sh```.
- This will create a subdirectory in *newFitData/human_SMARCAL1* with the same name as your project name, and populate it with as many coordinate files as you decided above in (iii).

*Note you will need to download the whole directory c++Molecule for this to run smoothly.*
