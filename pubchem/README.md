This collection consists of small, drug-like molecules from PubChem.  They were selected as follows.

First, SDF files were downloaded for all substances whose source was either BindingDB or ChemIDplus.  They were placed in the sources directory.  Due to the large size of the files, they have not been checked in.

processSDFFiles.py parses all files in that directory, applies a series of filters to select candidate molecules, and writes their IDs and SMILES strings to a corresponding text file in the same directory.

sortMolecules.py reads in the text files, combines them into a single list of molecules, and sorts them so that every molecule is maximally different from all the ones that have come before it.  It write the sorted molecules to sorted.txt.

createPubchem.py reads a subset of molecules from the list, generates 50 conformations for each one, and writes them to a HDF5 file.

The conformations for the first six sets (molecules 1-15000) were generated using OpenMM 7.6, RDKit 2020.09.3, OpenFF Toolkit 0.10.0, and MDTraj 1.9.5.

The conformations for the next four sets (molecules 15001-25000) were generated with OpenMM 8.0, RDKit 2023.03.1, OpenFF Toolkit 0.14.0, and MDTraj 1.9.7.

The conformations for molecules containing Boron and Silicon were generated with RDKit 2023.03.2, OpenFF Toolkit 0.14.0, and xtb 6.5.1.