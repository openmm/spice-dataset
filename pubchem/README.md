This collection consists of small, drug-like molecules from PubChem.  They were selected as follows.

First, SDF files were downloaded for all substances whose source was either BindingDB or ChemIDplus.  They were placed in the sources directory.  Due to the large size of the files, they have not been checked in.

processSDFFiles.py parses all files in that directory, applies a series of filters to select candidate molecules, and writes their IDs and SMILES strings to a corresponding text file in the same directory.

sortMolecules.py reads in the text files, combines them into a single list of molecules, and sorts them so that every molecule is maximally different from all the ones that have come before it.  It write the sorted molecules to sorted.txt.

createPubchem.py reads a subset of molecules from the list, generates 50 conformations for each one, and writes them to a HDF5 file.

The conformations were generated using OpenMM 7.6, RDKit 2020.09.3, OpenFF Toolkit 0.10.0, and MDTraj 1.9.5.