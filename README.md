# SPICE: A Dataset for Training Machine Learning Potentials

This repository contains scripts and data files used in the creation of the SPICE dataset.  It does not contain the
dataset itself.  That is available from Zenodo:

[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.10975225.svg)](https://doi.org/10.5281/zenodo.10975225)

SPICE (Small-Molecule/Protein Interaction Chemical Energies) is a collection of quantum mechanical data for
training potential functions.  The emphasis is particularly on simulating drug-like small molecules interacting
with proteins.  It is designed to achieve the following goals.

- **Cover a wide range of chemical space**.  It includes 17 elements (H, Li, B, C, N, O, F, Na, Mg, Si, P, S, Cl, K, Ca, Br, I)
  and a wide range of chemical groups.  It includes charged and polar molecules as well as neutral ones.  It is
  designed to sample a wide range of both covalent and non-covalent interactions.
- **Cover a wide range of conformations**.  It includes both low and high energy conformations.  It is
  designed to sample all regions of configuration space that are likely to be encountered in typical simulations.
- **Include forces as well as energies**.  Many datasets include only energies, not forces.  SPICE includes
  forces as well, which enormously increases the information content of the dataset.
- **Include a variety of other information**.  SPICE includes a variety of other QM results for each sample,
  such as bond orders, partial charges, and atomic multipoles.  It follows the principle that any result that
  is easy to compute and potentially useful should be made available.
- **Use an accuate level of theory**.  Models can never be more accurate than the data they are trained on.
  SPICE computations are done at the ωB97M-D3BJ/def2-TZVPPD level of theory.
- **Be a dynamic, growing dataset**.  The dataset should grow with time as new data is generated.  This will
  allow models trained on it both to improve in accuracy and to expand the range of chemical space they cover.
  Versioned releases will be created periodically to allow for reproducibility.
- **Be freely available under a non-restrictive licence**.  All data in the SPICE dataset may be used under the
  public domain equivalent [CC0 license](https://creativecommons.org/share-your-work/public-domain/cc0/).

SPICE is made up of a collection of subsets.  Each one is designed to provide a particular type of information.
The subsets in the current version (2.0) include the following.

- **Dipeptides**.  These provide comprehensive sampling of the covalent interactions found in proteins.
- **Solvated amino acids**.  These provide sampling of protein-water and water-water interactions.
- **PubChem molecules**.  These sample a very wide variety of drug-like small molecules.
- **Solvated PubChem molecules**.  These provide sampling of ligand-water interactions.
- **Monomer and dimer structures from [DES370K](https://www.nature.com/articles/s41597-021-00833-x)**.
  These provide sampling of a wide variety of non-covalent interactions.
- **Amino acid, ligand pairs**.  These provide sampling of nonbonded protein-ligand interactions.
- **Ion pairs**.  These provide further sampling of Coulomb interactions over a range of distances.
- **Water clusters**.  These provide additional sampling of water-water interactions.

This table summarizes the content of each subset: the number of molecules/clusters it contains, the total number of
conformations, the range of sizes spanned by the molecules/clusters, and the list of elements that appear in the subset.

|Subset|Molecules/Clusters|Conformations|Atoms|Elements|
|------|------------------|-------------|-----|--------|
|Dipeptides|677|33,850|26–60|H, C, N, O, S|
|Solvated Amino Acids|26|1300|79–96|H, C, N, O, S|
|DES370K Dimers|3490|345,676|2–34|H, Li, C, N, O, F, Na, Mg, P, S, Cl, K, Ca, Br, I|
|DES370K Monomers|374|18,700|3–22|H, C, N, O, F, P, S, Cl, Br, I|
|PubChem|28,039|1,398,566|3–50|H, B, C, N, O, F, Si, P, S, Cl, Br, I|
|Solvated PubChem|1397|13,934|63–110|H, C, N, O, F, P, S, Cl, Br, I|
|Amino Acid Ligand Pairs|79,967|194,174|24–72|H, C, N, O, F, P, S, Cl, Br, I|
|Ion Pairs|28|1426|2|Li, F, Na, Cl, K, Br, I|
|Water Clusters|1|1000|90|H, O|
|Total|113,999|2,008,628|2–110|H, Li, B, C, N, O, F, Na, Mg, Si, P, S, Cl, K, Ca, Br, I|

## Citing The Dataset

Please cite these manuscripts for papers that use the SPICE dataset:

Peter Eastman, Pavan Kumar Behara, David L. Dotson, Raimondas Galvelis, John E. Herr, Josh T. Horton, Yuezhi Mao,
John D. Chodera, Benjamin P. Pritchard, Yuanqing Wang, Gianni De Fabritiis, and Thomas E. Markland.  "SPICE, A Dataset
of Drug-like Molecules and Peptides for Training Machine Learning Potentials."  Scientific Data 10, 11 (2023).
https://doi.org/10.1038/s41597-022-01882-6

Peter Eastman, Benjamin P. Pritchard, John D. Chodera, Thomas E. Markland.  "Nutmeg and SPICE: Models and Data for
Biomolecular Machine Learning."  J. Chem. Theory Comput. 20, 19, 8583-8593 (2024).  https://doi.org/10.1021/acs.jctc.4c00794

To cite a particular version of the dataset, cite the Zenodo DOI found on the Releases page and shown above for the
most recent version.

## Generating New Data

All calculations in the SPICE dataset are computed with [Psi4](https://github.com/psi4/psi4).  If you want to generate
new data that can be combined with SPICE, it is important to use the same level of theory and the same program with
identical settings.  Even when two programs use the same level of theory, there usually are enough differences in how
they do calculations that energies they produce cannot be directly compared to each other.  A
[sample input file](sample.dat) for Psi4 is provided.  It  shows the exact settings to use to produce new data that can
be correctly combined with SPICE.