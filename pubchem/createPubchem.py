import openmm
import openmm.app as app
import openmm.unit as unit
from openff.toolkit.typing.engines.smirnoff import ForceField
from openff.toolkit.topology import Molecule, Topology
from concurrent.futures import ThreadPoolExecutor
import numpy as np
import h5py
import sys

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

def createConformations(outputfile, forcefield, smiles, sid):
    """Generate the conformations for a molecule and save them to disk."""
    print(f'Generating {index}: {smiles}')
    try:
        mol = Molecule.from_smiles(smiles, allow_undefined_stereo=True)
        fftop = Topology()
        fftop.add_molecule(mol)
        mmtop = fftop.to_openmm()
        system = forcefield.create_openmm_system(fftop)
    except:
        print('  failed to parametrize')
        return

    # Generate 10 diverse starting points.  Run MD from each one to generate a total
    # of 100 high energy conformations.

    mol.generate_conformers(n_conformers=10, rms_cutoff=0*unit.nanometers)
    assert len(mol.conformers) == 10
    
    def simulate(pos):
        integrator = openmm.LangevinMiddleIntegrator(500*unit.kelvin, 1/unit.picosecond, 0.001*unit.picosecond)
        simulation = app.Simulation(mmtop, system, integrator, openmm.Platform.getPlatformByName('Reference'))
        simulation.context.setPositions(pos)
        simulation.minimizeEnergy()
        simulation.context.setVelocitiesToTemperature(500*unit.kelvin)
        states = []
        for i in range(10):
            simulation.step(10000)
            state = simulation.context.getState(getPositions=True, getEnergy=True)
            if state.getPotentialEnergy() < 1e4*unit.kilojoules_per_mole:
                states.append(state)
        return states

    futures = []
    with ThreadPoolExecutor() as executor:
        for pos in mol.conformers:
            futures.append(executor.submit(simulate, pos))
    states = []
    for future in futures:
        states += future.result()

    # Select 25 that are most different from each other.

    if len(states) < 25:
        print('  failed to generate states')
        return
    states = filterByRMSD(states, mmtop)

    # Create a nearby, lower energy conformation from each one.

    integrator = openmm.LangevinMiddleIntegrator(100*unit.kelvin, 1/unit.picosecond, 0.001*unit.picosecond)
    simulation = app.Simulation(mmtop, system, integrator, openmm.Platform.getPlatformByName('Reference'))
    for state in states[:]:
        simulation.context.setState(state)
        simulation.minimizeEnergy(maxIterations=5)
        simulation.context.setVelocitiesToTemperature(100*unit.kelvin)
        simulation.step(1000)
        states.append(simulation.context.getState(getPositions=True))
    saveToFile(outputfile, mol, states, sid)

def saveToFile(outputfile, mol, states, name):
    """Save a molecule and its conformations to a HDF5 file."""
    try:
        mol._conformers = None
        for state in states:
            mol.add_conformer(state.getPositions(asNumpy=True))
        mol = mol.canonical_order_atoms()
        smiles = mol.to_smiles(isomeric=True, explicit_hydrogens=True, mapped=True)
    except:
        print('  exception generating canonical SMILES')
        return
    conformations = [c.value_in_unit(unit.nanometers) for c in mol.conformers]
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

first = int(sys.argv[1])
last = int(sys.argv[2])
print(f'Creating conformations for molecules {first} to {last}')
forcefield = ForceField('openff_unconstrained-2.0.0.offxml')
outputfile = h5py.File(f'pubchem-{first}-{last}.hdf5', 'w')
index = 0
for line in open('sorted.txt'):
    index += 1
    if index < first:
        continue
    if index > last:
        break
    sid, smiles = line.split()
    createConformations(outputfile, forcefield, smiles, sid)
