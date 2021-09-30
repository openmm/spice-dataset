import openmm.unit as unit
import numpy as np
import dask.dataframe as dd
import h5py
import os
from collections import defaultdict
from rdkit import Chem
from openff.toolkit.topology import Molecule

# Load the SDF files for the monomers.

print('Loading monomers')
moleculeForSmiles = {}
for filename in os.listdir('SDFS'):
    if not filename.endswith('.sdf'):
        continue
    smiles = filename[:-4]
    supp = Chem.SDMolSupplier(f'SDFS/{filename}', sanitize=False, removeHs=False)
    moleculeForSmiles[smiles] = list(supp)[0]

# Find all rows for each unique pair of monomers.

print('Collating samples')
skip = ['[Ar]', '[He]', '[Kr]', '[Ne]', '[Xe]']
coordsForDimer = defaultdict(list)
df = dd.read_csv(f'Donchev et al DES370K.csv')
for i, row in enumerate(df.itertuples()):
    smiles1 = row.smiles0
    smiles2 = row.smiles1
    if smiles1 in skip or smiles2 in skip:
        continue
    rdmol1 = moleculeForSmiles[smiles1]
    rdmol2 = moleculeForSmiles[smiles2]
    symbols = [a.GetSymbol() for a in rdmol1.GetAtoms()] + [a.GetSymbol() for a in rdmol2.GetAtoms()]
    expectedElements = ' '.join(symbols)
    assert str(row.elements) == expectedElements
    coords = np.array([float(f) for f in row.xyz.split()]).reshape(len(symbols), 3)
    coordsForDimer[(smiles1, smiles2)].append(coords)

# Process the dimers and create the output file.

print('Creating file')
outputfile = h5py.File('des370K.hdf5', 'w')
for dimer in coordsForDimer:
    smiles1, smiles2 = dimer
    mol1 = Molecule.from_rdkit(moleculeForSmiles[smiles1], allow_undefined_stereo=True, hydrogens_are_explicit=True)
    mol2 = Molecule.from_rdkit(moleculeForSmiles[smiles2], allow_undefined_stereo=True, hydrogens_are_explicit=True)

    # Put the atoms and coordinates into canonical order.

    mol1._conformers = None
    mol2._conformers = None
    for coords in coordsForDimer[dimer]:
        mol1.add_conformer(coords[:mol1.n_atoms]*unit.angstrom)
        mol2.add_conformer(coords[mol1.n_atoms:]*unit.angstrom)
    mol1 = mol1.canonical_order_atoms()
    mol2 = mol2.canonical_order_atoms()
    conformations = [np.concatenate([c1.value_in_unit(unit.nanometers), c2.value_in_unit(unit.nanometers)], axis=0)
                     for c1, c2 in zip(mol1.conformers, mol2.conformers)]
    conformations = np.array(conformations)
    assert conformations.shape == (len(coordsForDimer[dimer]), mol1.n_atoms+mol2.n_atoms, 3)

    # Compute the SMILES string.

    s1 = mol1.to_smiles(isomeric=True, explicit_hydrogens=True, mapped=True)
    mapping = dict((i, i+mol1.n_atoms+1) for i in range(mol2.n_atoms))
    mol2.properties["atom_map"] = mapping
    s2 = mol2.to_smiles(isomeric=True, explicit_hydrogens=True, mapped=True)
    smiles = f'{s1}.{s2}'

    # Write the data to the output file.

    group = outputfile.create_group(f'{smiles1} {smiles2}')
    group.create_dataset('smiles', data=[smiles], dtype=h5py.string_dtype())
    ds = group.create_dataset('conformations', data=conformations, dtype=np.float32)
    ds.attrs['units'] = 'nanometers'
