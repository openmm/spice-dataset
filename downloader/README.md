# Data Downloader

The full dataset is available on [QCArchive](https://qcarchive.molssi.org/).  It can be downloaded using the [QCPortal](https://docs.qcarchive.molssi.org/projects/QCPortal/en/stable/) library.  Although QCPortal is flexible, it is not very convenient for most purposes.  This directory contains a script to download exactly the information you want and store it all in a single HDF5 file.

To run it, you first need to install a few libraries that it requires.  Use this command.

    conda install -c conda-forge qcportal pyyaml h5py rdkit numpy

Next edit the `config.yaml` file.  Select the data subsets to download (most often you will want all of them) and data values to include in the output file (`'DFT TOTAL ENERGY'` and `'DFT TOTAL GRADIENT'` are the most common ones; others are only needed for special purposes).  Comment out any subset or value you do not want by placing a `#` at the start of the line.

Finally run `python downloader.py` to download the data.  This can take several hours.  The data is stored in a HDF5 file called `SPICE.hdf5`.  It is structured as follows.

- There is one top level group for each unique molecule or cluster.  The name of each group is either a PubChem Substance ID (for PubChem molecules), an amino acid sequence (for dipeptides and solvated amino acids), or a SMILES string (for everything else).
- Each group contains the following datasets.  `N` is the number of atoms in the molecule and `M` is the number of conformations.
  - `subset`: The name of the data subset the molecule is from.
  - `smiles`: The canonical SMILES string for the molecule.  It includes explicit hydrogens and atom indices.
  - `atomic_numbers`: Array of length `N` containing the atomic number of every atom.  They are ordered following the indices in the SMILES string.
  - `conformations`: Array of shape `(M, N, 3)` containing the atomic coordinates for every conformation.
  - `formation_energy`: Array of length `M` containing the total energy of each conformation, minus the reference energies of the individual atoms when infinitely separated.  This is the most useful energy for most purposes, since it contains all energy components that vary with atom positions but removes the large constant part corresponding to the internal energies of individual atoms.
  - Other values as specified in the configuration file.  For example, `dft_total_energy` is an array of length `M` containing the energy of each conformation.  `dft_total_gradient` is an array of shape `(M, N, 3)` containing the gradient of the energy with respect to the atomic coordinates.
- All values are in atomic units.  Distances are in bohr and energies in hartree.