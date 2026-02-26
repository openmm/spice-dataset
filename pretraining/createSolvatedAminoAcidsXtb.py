import openmm
import openmm.app as app
import openmm.unit as unit
import numpy as np
import h5py
from utils import compute_xtb_for_conformers, convert_to_openff, save_to_file
from collections import namedtuple

def createModel(forcefield, res, var):
    """Create a capped amino acid of the specified type and build a water box around it."""
    import pdbfixer
    fixer = pdbfixer.PDBFixer(filename='../solvated-amino-acids/ace_ala_nme.pdb')
    fixer.missingResidues = {}
    fixer.applyMutations([f'ALA-2-{res}'], 'A')
    fixer.findMissingAtoms()
    fixer.addMissingAtoms()
    modeller = app.Modeller(fixer.topology, fixer.positions)
    modeller.addHydrogens(forcefield=forcefield, variants=[None, var, None])
    modeller.addSolvent(forcefield, boxSize=(2.2, 2.2, 2.2)*unit.nanometer)
    return modeller.topology, modeller.positions

def createConformations(outputfile, forcefield, topology, positions, name, charges):
    """Generate the conformations for a molecule and save them to disk."""
    print(f'Generating {name}')
    system = forcefield.createSystem(topology, nonbondedMethod=app.CutoffPeriodic, constraints=None, rigidWater=False)
    system.addForce(openmm.MonteCarloBarostat(1*unit.bar, 300*unit.kelvin))
    integrator = openmm.LangevinMiddleIntegrator(300*unit.kelvin, 1/unit.picosecond, 0.001*unit.picosecond)
    simulation = app.Simulation(topology, system, integrator)
    simulation.context.setPositions(positions)

    # Equilibrate for 100 ps.

    simulation.minimizeEnergy()
    simulation.context.setVelocitiesToTemperature(300*unit.kelvin)
    simulation.step(100000)

    # Simulate for 5 ns, saving a state every 10 ps.

    import mdtraj
    atoms = list(topology.atoms())
    soluteAtomIndex = [atom.index for atom in atoms if atom.residue.chain.index == 0]
    waterAtomIndex = [atom.index for atom in atoms if atom.residue.name == 'HOH' and atom.name == 'O']
    from openmm.app.internal import compiled
    conformations = []
    for i in range(509):
        simulation.step(10000)
        state = simulation.context.getState(getPositions=True, getEnergy=True)
        assert state.getPotentialEnergy() < 0*unit.kilojoules_per_mole

        # Identify the 20 water molecules closest to the solute.

        periodicDistance = compiled.periodicDistance(state.getPeriodicBoxVectors().value_in_unit(unit.nanometer))
        pos = state.getPositions().value_in_unit(unit.nanometer)
        distances = [min(periodicDistance(pos[i], pos[j]) for j in soluteAtomIndex) for i in waterAtomIndex]
        sortedIndices = np.argsort(distances)
        waterToKeep = set(atoms[waterAtomIndex[i]].residue.index for i in sortedIndices[:20])

        # Remove everything else from the Topology and coordinates.

        modeller = app.Modeller(topology, state.getPositions())
        toDelete = [res for res in topology.residues() if res.chain.index > 0 and res.index not in waterToKeep]
        modeller.delete(toDelete)

        # Center the solute with the water around it.

        mdtop = mdtraj.Topology.from_openmm(modeller.topology)
        xyz = np.array([modeller.positions.value_in_unit(unit.nanometer)])
        traj = mdtraj.Trajectory(xyz, mdtop)
        traj.unitcell_vectors = np.array([state.getPeriodicBoxVectors().value_in_unit(unit.nanometer)])
        traj = traj.image_molecules()
        conformations.append(traj.xyz[0]*unit.nanometer)

    # Create an OpenFF molecule.

    system = forcefield.createSystem(modeller.topology, constraints=None, rigidWater=False)
    integrator = openmm.LangevinMiddleIntegrator(300*unit.kelvin, 1/unit.picosecond, 0.001*unit.picosecond)
    simulation = app.Simulation(modeller.topology, system, integrator)
    simulation.context.setPositions(conformations[-1])
    mol = convert_to_openff(simulation, charges)

    # Add the conformations to the Molecule and put the atoms in canonical order.

    mol._conformers = None
    for conf in conformations:
        mol.add_conformer(conf)
    mol = mol.canonical_order_atoms()

    # Compute energies and forces with XTB.

    positions, energies, formation_energies, grads = compute_xtb_for_conformers(mol)
    return positions, energies, formation_energies, grads, mol

# Define the residue variants we will include.

Residue = namedtuple('Residue', ['name', 'variant', 'charges'])

residues = [
    Residue('ALA', None, {}),
    Residue('ASN', None, {}),
    Residue('CYS', 'CYS', {}),
    Residue('CYS', 'CYX', {'SG':[-1,True]}),
    Residue('GLU', 'GLH', {}),
    Residue('GLU', 'GLU', {'OE1':[-1,False], 'OE2':[-1,False]}),
    Residue('HIS', 'HID', {}),
    Residue('HIS', 'HIE', {}),
    Residue('HIS', 'HIP', {'ND1':[1,True]}),
    Residue('LEU', None, {}),
    Residue('MET', None, {}),
    Residue('PRO', None, {}),
    Residue('THR', None, {}),
    Residue('TYR', None, {}),
    Residue('ARG', None, {'NH1':[1,False], 'NH2':[1,False]}),
    Residue('ASP', 'ASH', {}),
    Residue('ASP', 'ASP', {'OD1':[-1,False], 'OD2':[-1,False]}),
    Residue('GLN', None, {}),
    Residue('GLY', None, {}),
    Residue('ILE', None, {}),
    Residue('LYS', 'LYN', {}),
    Residue('LYS', 'LYS', {'NZ':[1,True]}),
    Residue('PHE', None, {}),
    Residue('SER', None, {}),
    Residue('TRP', None, {}),
    Residue('VAL', None, {})
]

# Create the conformations for the molecules.

forcefield = app.ForceField('amber14-all.xml', 'amber14/tip3pfb.xml')
outputfile = h5py.File('solvated-amino-acids.hdf5', 'w')
for res in residues:
    name = res.name if res.variant is None else res.variant
    topology, positions = createModel(forcefield, res.name, res.variant)
    charges = {}
    for atom in topology.atoms():
        if atom.residue.index == 1 and atom.name in res.charges:
            charges[atom.index] = res.charges[atom.name]
    positions, energies, formation_energies, grads, mol = createConformations(outputfile, forcefield, topology, positions, name, charges)
    save_to_file(outputfile, mol, positions, energies, formation_energies, grads, name, 'Solvated Amino Acids')
