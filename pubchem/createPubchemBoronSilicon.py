import ase
import ase.md
import ase.md.velocitydistribution
import ase.optimize
from xtb.ase.calculator import XTB
from openff.toolkit.topology import Molecule
from openff.units import unit as ffunit
from concurrent.futures import ProcessPoolExecutor
import numpy as np
import h5py
import os
import logging

logging.getLogger().setLevel(logging.ERROR)

def filterByRMSD(states, mol):
    """From a set of States, return 25 whose RMSDs to each other are maximally different."""
    import mdtraj
    xyz = np.array(states)
    top = mdtraj.Topology()
    chain = top.add_chain()
    res = top.add_residue('mol', chain)
    for atom in mol.atoms:
        top.add_atom(atom.name, mdtraj.element.Element.getByAtomicNumber(atom.atomic_number), res)
    traj = mdtraj.Trajectory(xyz, top)
    traj.center_coordinates()
    finalStates = set([0])
    minRmsd = mdtraj.rmsd(traj, traj, 0, precentered=True)
    for i in range(24):
        best = np.argmax(minRmsd)
        minRmsd = np.minimum(minRmsd, mdtraj.rmsd(traj, traj, best, precentered=True))
        finalStates.add(best)
    return [states[i] for i in finalStates]

def simulate(pos, numbers, charges):
    atoms = ase.Atoms(positions=pos, numbers=numbers, charges=charges)
    atoms.calc = XTB(method="GFN-FF")
    ase.optimize.LBFGS(atoms, logfile=os.devnull).run(0.001, 20)
    ase.md.velocitydistribution.MaxwellBoltzmannDistribution(atoms, temperature_K=300)
    dyn = ase.md.langevin.Langevin(atoms, 1*ase.units.fs, temperature_K=300, friction=1e-3)
    states = []
    for i in range(10):
        dyn.run(10000)
        states.append((atoms.get_positions()-np.mean(atoms.get_positions(), axis=0)))
    return states

def createConformations(smiles, sid, mol):
    """Generate the conformations for a molecule."""
    print(f'Generating {smiles}')
    numbers = [a.atomic_number for a in mol.atoms]
    charges = [a.formal_charge.m for a in mol.atoms]

    # Generate 10 diverse starting points.  Run MD from each one to generate a total
    # of 100 high energy conformations.

    mol.generate_conformers(n_conformers=10, rms_cutoff=0*ffunit.nanometers)
    assert len(mol.conformers) == 10
    states = []
    for pos in mol.conformers:
        states += simulate(pos.m_as(ffunit.angstroms), numbers, charges)

    # Select 25 that are most different from each other.

    if len(states) < 25:
        print('  failed to generate states')
        return
    states = filterByRMSD(states, mol)

    # Create a nearby, lower energy conformation from each one.

    for state in states[:]:
        atoms = ase.Atoms(positions=state, numbers=numbers, charges=charges)
        atoms.calc = XTB(method="GFN-FF")
        ase.optimize.LBFGS(atoms, logfile=os.devnull).run(0.001, 5)
        ase.md.velocitydistribution.MaxwellBoltzmannDistribution(atoms, temperature_K=100)
        dyn = ase.md.langevin.Langevin(atoms, 1*ase.units.fs, temperature_K=100, friction=1e-3)
        dyn.run(1000)
        states.append(atoms.get_positions())
    return mol, states, sid

def saveToFile(outputfile, mol, states, name):
    """Save a molecule and its conformations to a HDF5 file."""
    try:
        mol._conformers = None
        for state in states:
            mol.add_conformer(state*ffunit.angstroms)
        mol = mol.canonical_order_atoms()
        smiles = mol.to_smiles(isomeric=True, explicit_hydrogens=True, mapped=True)
    except:
        print('  exception generating canonical SMILES')
        return
    conformations = [c.m_as(ffunit.nanometers) for c in mol.conformers]
    conformations = [c-np.average(c, axis=0) for c in conformations]
    group = outputfile.create_group(name)
    group.create_dataset('smiles', data=[smiles], dtype=h5py.string_dtype())
    ds = group.create_dataset('conformations', data=np.array(conformations), dtype=np.float32)
    ds.attrs['units'] = 'nanometers'

    # As a sanity check, make sure the SMILES string doesn't have any radicals.

    from rdkit import Chem
    rdmol = Chem.MolFromSmiles(smiles)
    assert all(atom.GetNumRadicalElectrons() == 0 for atom in rdmol.GetAtoms())

# Create the molecules.

outputfile = h5py.File(f'pubchem-boron-silicon.hdf5', 'w')
futures = []
with ProcessPoolExecutor() as executor:
    for filename in ['BindingDB1.txt', 'BindingDB2.txt', 'BindingDB3.txt', 'ChemIDplus.txt']:
        for line in open('sources/'+filename):
            sid, smiles = line.split()
            mol = Molecule.from_smiles(smiles, allow_undefined_stereo=True)
            if any(a.symbol in ('B', 'Si') for a in mol.atoms):
                futures.append(executor.submit(createConformations, smiles, sid, mol))

# Save the results to the output file.

for future in futures:
    mol, states, sid = future.result()
    saveToFile(outputfile, mol, states, sid)