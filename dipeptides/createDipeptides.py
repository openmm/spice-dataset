import openmm
import openmm.app as app
import openmm.unit as unit
import numpy as np
import h5py
from collections import namedtuple

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

def convertToOpenFF(simulation, charges):
    """Create an OpenForceField Molecule object from an OpenMM simulation."""
    # Use Open Babel to identify bond orders and stereochemistry.

    from openbabel import openbabel
    obConversion = openbabel.OBConversion()
    obConversion.SetInFormat("pdb")
    from io import StringIO
    io = StringIO()
    positions = simulation.context.getState(getPositions=True).getPositions()
    app.PDBFile.writeFile(simulation.topology, positions, io)
    obmol = openbabel.OBMol()
    obConversion.ReadString(obmol, io.getvalue())
    facade = openbabel.OBStereoFacade(obmol)

    # Build the molecule.

    from openff.toolkit.topology import Molecule
    mol = Molecule()
    for i, atom in enumerate(simulation.topology.atoms()):
        obatom = obmol.GetAtom(i+1)
        charge = 0
        if atom.index in charges:
            charge = charges[atom.index]
        stereo = None
        if facade.HasTetrahedralStereo(obatom.GetId()):
            stereo = ('S', 'R')[facade.GetTetrahedralStereo(obatom.GetId()).GetConfig().winding]
        mol.add_atom(atom.element.atomic_number, charge, obatom.IsAromatic(), stereo, atom.name)
    for i in range(obmol.NumBonds()):
        bond = obmol.GetBond(i)
        mol.add_bond(bond.GetBeginAtomIdx()-1, bond.GetEndAtomIdx()-1, bond.GetBondOrder(), bond.IsAromatic())

    # Verify that the total charge is correct.

    nb = [f for f in simulation.context.getSystem().getForces() if isinstance(f, openmm.NonbondedForce)][0]
    charge = round(sum(nb.getParticleParameters(i)[0]._value for i in range(nb.getNumParticles())))
    assert charge == mol.total_charge.value_in_unit(unit.elementary_charge)
    return mol

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

def createConformations(outputfile, forcefield, topology, positions, name, charges):
    """Generate the conformations for a molecule and save them to disk."""
    print(f'Generating {name}')
    system = forcefield.createSystem(topology, constraints=None, rigidWater=False)
    integrator = openmm.LangevinMiddleIntegrator(500*unit.kelvin, 1/unit.picosecond, 0.001*unit.picosecond)
    simulation = app.Simulation(topology, system, integrator, openmm.Platform.getPlatformByName('Reference'))
    simulation.context.setPositions(positions)
    simulation.minimizeEnergy()

    # Generate 10 diverse starting points.  Run MD from each one to generate a total
    # of 100 high energy conformations.

    mol = convertToOpenFF(simulation, charges)
    mol.generate_conformers(n_conformers=10, rms_cutoff=0*unit.nanometers)
    assert len(mol.conformers) == 10
    states = []
    for pos in mol.conformers:
        simulation.context.setPositions(pos)
        simulation.minimizeEnergy()
        simulation.context.setVelocitiesToTemperature(500*unit.kelvin)
        for i in range(10):
            simulation.step(10000)
            states.append(simulation.context.getState(getPositions=True))

    # Select 25 that are most different from each other.

    states = filterByRMSD(states, topology)

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

# Define the residue variants we will include.

Residue = namedtuple('Residue', ['name', 'variant', 'charges'])

residues = [
    Residue('ALA', None, {}),
    Residue('ASN', None, {}),
    Residue('CYS', 'CYS', {}),
    Residue('CYS', 'CYX', {'SG':-1}),
    Residue('GLU', 'GLH', {}),
    Residue('GLU', 'GLU', {'OE1':-1}),
    Residue('HIS', 'HID', {}),
    Residue('HIS', 'HIE', {}),
    Residue('HIS', 'HIP', {'ND1':1}),
    Residue('LEU', None, {}),
    Residue('MET', None, {}),
    Residue('PRO', None, {}),
    Residue('THR', None, {}),
    Residue('TYR', None, {}),
    Residue('ARG', None, {'NH1':1}),
    Residue('ASP', 'ASH', {}),
    Residue('ASP', 'ASP', {'OD1':-1}),
    Residue('GLN', None, {}),
    Residue('GLY', None, {}),
    Residue('ILE', None, {}),
    Residue('LYS', 'LYN', {}),
    Residue('LYS', 'LYS', {'NZ':1}),
    Residue('PHE', None, {}),
    Residue('SER', None, {}),
    Residue('TRP', None, {}),
    Residue('VAL', None, {})
]

# Create the molecules.

forcefield = app.ForceField('amber14-all.xml')
outputfile = h5py.File('dipeptides.hdf5', 'w')
pdb = app.PDBFile('disulfide.pdb')
createConformations(outputfile, forcefield, pdb.topology, pdb.positions, 'Disulfide', {})
for res1 in residues:
    name1 = res1.name if res1.variant is None else res1.variant
    for res2 in residues:
        name2 = res2.name if res2.variant is None else res2.variant
        topology, positions = createDipeptide(forcefield, res1.name, res2.name, res1.variant, res2.variant)
        charges = {}
        for atom in topology.atoms():
            if atom.residue.index == 1 and atom.name in res1.charges:
                charges[atom.index] = res1.charges[atom.name]
            elif atom.residue.index == 2 and atom.name in res2.charges:
                charges[atom.index] = res2.charges[atom.name]
        createConformations(outputfile, forcefield, topology, positions, f'{name1}-{name2}', charges)
