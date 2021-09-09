import openmm
import openmm.app as app
import openmm.unit as unit
import numpy as np
import h5py

def createDipeptide(forcefield, res1, res2, var1, var2):
    """Use PDBFixer to create a dipeptide with specified residues and variants."""
    import pdbfixer
    fixer = pdbfixer.PDBFixer(filename='ala_ala.pdb')
    fixer.missingResidues = {}
    fixer.applyMutations([f'ALA-2-{res1}', f'ALA-3-{res2}'], 'A')
    fixer.findMissingAtoms()
    fixer.addMissingAtoms()
    modeller = app.Modeller(fixer.topology, fixer.positions)
    modeller.addHydrogens(forcefield=forcefield, variants=[None, var1, var2, None])
    return modeller.topology, modeller.positions

def convertToRdkit(simulation):
    """Create a RDKit Mol object from an OpenMM simulation."""
    from rdkit import Chem
    from io import StringIO
    io = StringIO()
    positions = simulation.context.getState(getPositions=True).getPositions()
    app.PDBFile.writeFile(simulation.topology, positions, io)
    mol = Chem.rdmolfiles.MolFromPDBBlock(io.getvalue(), sanitize=False, removeHs=False)
    mol.RemoveAllConformers()
    return mol

def findRdkitConformers(simulation):
    """Use RDKit to generate 10 diverse conformations for a molecule."""
    from rdkit.Chem import AllChem
    mol = convertToRdkit(simulation)
    AllChem.EmbedMultipleConfs(mol, numConfs=10)
    conformers = [mol.GetConformer(i).GetPositions()*unit.angstrom for i in range(mol.GetNumConformers())]
    return conformers

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

def createConformations(outputfile, forcefield, topology, positions, name):
    """Generate the conformations for a molecule and save them to disk."""
    print(f'Generating {name}')
    system = forcefield.createSystem(topology, constraints=None, rigidWater=False)
    integrator = openmm.LangevinMiddleIntegrator(500*unit.kelvin, 1/unit.picosecond, 0.001*unit.picosecond)
    simulation = app.Simulation(topology, system, integrator, openmm.Platform.getPlatformByName('Reference'))
    simulation.context.setPositions(positions)
    simulation.minimizeEnergy()

    # Use RDKit to generate 10 diverse starting points.  Run MD from each one
    # to generate a total of 100 high energy conformations.

    startPositions = findRdkitConformers(simulation)
    states = []
    for pos in startPositions:
        simulation.context.setPositions(pos)
        simulation.minimizeEnergy()
        simulation.context.setVelocitiesToTemperature(500*unit.kelvin)
        for i in range(10):
            simulation.step(10000)
            states.append(simulation.context.getState(getPositions=True))

    # Select 25 that are most different from each other.

    states = filterByRMSD(states, topology)

    # Do a few iterations of energy minimization to create a nearby, lower
    # energy conformation from each one.

    for state in states[:]:
        simulation.context.setState(state)
        simulation.minimizeEnergy(maxIterations=3)
        states.append(simulation.context.getState(getPositions=True))
    saveToFile(outputfile, simulation, states, name)

def saveToFile(outputfile, simulation, states, name):
    """Save a molecule and its conformations to a HDF5 file."""
    from openff.toolkit.topology import Molecule
    rdmol = convertToRdkit(simulation)
    mol = Molecule.from_rdkit(rdmol, allow_undefined_stereo=True, hydrogens_are_explicit=True)
    for state in states:
        mol.add_conformer(state.getPositions(asNumpy=True))
    mol = mol.canonical_order_atoms()
    smiles = mol.to_smiles(isomeric=True, explicit_hydrogens=True, mapped=True)
    conformations = [c.value_in_unit(unit.nanometers) for c in mol.conformers]
    group = outputfile.create_group(name)
    group.create_dataset('smiles', data=[smiles], dtype=h5py.string_dtype())
    ds = group.create_dataset('conformations', data=np.array(conformations), dtype=np.float32)
    ds.attrs['units'] = 'nanometers'

# Create the molecules.

forcefield = app.ForceField('amber14-all.xml')
outputfile = h5py.File('dipeptides.hdf5', 'w')
pdb = app.PDBFile('disulfide.pdb')
createConformations(outputfile, forcefield, pdb.topology, pdb.positions, 'Disulfide')
residues = ['ALA', 'ASN', 'CYS', 'CYS', 'GLU', 'GLU', 'HIS', 'HIS', 'HIS', 'LEU', 'MET', 'PRO', 'THR', 'TYR', 'ARG', 'ASP', 'ASP', 'GLN', 'GLY', 'ILE', 'LYS', 'LYS', 'PHE', 'SER', 'TRP', 'VAL']
variants = [None,  None,  'CYS', 'CYX', 'GLH', 'GLU', 'HID', 'HIE', 'HIP', None,  None,  None,  None,  None,  None,  'ASH', 'ASP', None,  None,  None,  'LYN', 'LYS', None,  None,  None,  None]
for res1, var1 in zip(residues, variants):
    name1 = res1 if var1 is None else var1
    for res2, var2 in zip(residues, variants):
        name2 = res2 if var2 is None else var2
        topology, positions = createDipeptide(forcefield, res1, res2, var1, var2)
        createConformations(outputfile, forcefield, topology, positions, f'{name1}-{name2}')
