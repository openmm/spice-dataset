import openmm
import openmm.app as app
import openmm.unit as unit
import numpy as np
import h5py
import os
from rdkit import Chem
from openff.toolkit.typing.engines.smirnoff import ForceField
from openff.toolkit.topology import Molecule, Topology

def filterByRMSD(states, topology):
    """From a set of States, return 25 whose RMSDs to each other are maximally different."""
    import mdtraj
    xyz = np.array([s.getPositions().value_in_unit(unit.nanometer) for s in states])
    traj = mdtraj.Trajectory(xyz, mdtraj.Topology.from_openmm(topology))
    traj.center_coordinates()
    finalStates = set([0])
    minRmsd = mdtraj.rmsd(traj, traj, 0, precentered=True)
    for i in range(24):
        best = np.argmax(minRmsd)
        minRmsd = np.minimum(minRmsd, mdtraj.rmsd(traj, traj, best, precentered=True))
        finalStates.add(best)
    return [states[i] for i in finalStates]

def createConformations(outputfile, forcefield, mol, name):
    """Generate the conformations for a molecule and save them to disk."""
    print(f'Generating {name}')
    topology = Topology.from_molecules([mol])
    mmTopology = topology.to_openmm()
    system = forcefield.create_openmm_system(topology)
    system.addForce(openmm.CMMotionRemover())
    integrator = openmm.LangevinMiddleIntegrator(500*unit.kelvin, 1/unit.picosecond, 0.001*unit.picosecond)
    simulation = app.Simulation(mmTopology, system, integrator, openmm.Platform.getPlatformByName('Reference'))

    # Generate 10 diverse starting points.  Run MD from each one to generate a total
    # of 100 high energy conformations.

    mol.generate_conformers(n_conformers=10, rms_cutoff=0*unit.nanometers)
    assert len(mol.conformers) == 10
    states = []
    for pos in mol.conformers:
        simulation.context.setPositions(pos)
        simulation.minimizeEnergy()
        simulation.context.setVelocitiesToTemperature(500*unit.kelvin)
        for i in range(10):
            simulation.step(10000)
            state = simulation.context.getState(getPositions=True, getEnergy=True)
            if state.getPotentialEnergy() < 1e4*unit.kilojoules_per_mole:
                states.append(state)

    # Select 25 that are most different from each other.

    states = filterByRMSD(states, mmTopology)

    # Create a nearby, lower energy conformation from each one.

    integrator.setTemperature(100*unit.kelvin)
    for state in states[:]:
        simulation.context.setState(state)
        simulation.minimizeEnergy(maxIterations=5)
        simulation.context.setVelocitiesToTemperature(100*unit.kelvin)
        simulation.step(1000)
        states.append(simulation.context.getState(getPositions=True))
    saveToFile(outputfile, mol, states, name)

def saveToFile(outputfile, mol, states, name):
    """Save a molecule and its conformations to a HDF5 file."""
    mol._conformers = None
    for state in states:
        mol.add_conformer(state.getPositions(asNumpy=True))
    mol = mol.canonical_order_atoms()
    smiles = mol.to_smiles(isomeric=True, explicit_hydrogens=True, mapped=True)
    conformations = [c.value_in_unit(unit.nanometers) for c in mol.conformers]
    group = outputfile.create_group(name)
    group.create_dataset('smiles', data=[smiles], dtype=h5py.string_dtype())
    ds = group.create_dataset('conformations', data=np.array(conformations), dtype=np.float32)
    ds.attrs['units'] = 'nanometers'

    # As a sanity check, make sure the SMILES string doesn't have any radicals.

    rdmol = Chem.MolFromSmiles(smiles)
    assert all(atom.GetNumRadicalElectrons() == 0 for atom in rdmol.GetAtoms())

# Create the molecules.

forcefield = ForceField('openff_unconstrained-1.3.0.offxml')
outputfile = h5py.File('des-monomers.hdf5', 'w')
for filename in os.listdir('SDFS'):
    if not filename.endswith('.sdf'):
        continue
    smiles = filename[:-4]
    supp = Chem.SDMolSupplier(f'SDFS/{filename}', sanitize=False, removeHs=False)
    rdmol = list(supp)[0]
    if rdmol.GetNumAtoms() > 1:
        try:
            mol = Molecule.from_rdkit(rdmol, allow_undefined_stereo=True, hydrogens_are_explicit=True)
            createConformations(outputfile, forcefield, mol, smiles)
        except:
            print('  failed to parametrize')
