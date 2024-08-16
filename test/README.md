This is a test set to be used for evaluating models trained on the SPICE dataset.  It should not be
used for training.  It consists of 600 molecules that are not present in the training set, with 10
conformations for each molecule.

- 200 molecules from Ligand Expo with between 40 and 50 atoms.
- 200 molecules from Ligand Expo with between 70 and 80 atoms.
- 200 pentapeptides with up to 110 atoms.

The conformations were generated using OpenMM 8.1.1, RDKit 2023.09.6, OpenFF Toolkit 0.15.2, ASE 3.22.1,
and xtb 6.5.1.

A second test set consists of 200 ligand/amino-acid dimers with 10 conformations for each one.
They consist of molecules that are present in the training set, but not in the same combinations.