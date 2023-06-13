import openmm
import openmm.app as app
import openmm.unit as unit
import h5py
import os
from utils import compute_xtb_for_conformers, convert_to_openff, save_to_file
from openff.units import unit as ffunit
from collections import namedtuple
from concurrent.futures import ThreadPoolExecutor

def createDipeptide(forcefield, res1, res2, var1, var2):
    """Use PDBFixer to create a dipeptide with specified residues and variants."""
    import pdbfixer
    fixer = pdbfixer.PDBFixer(filename='../dipeptides/ala_ala.pdb')
    fixer.missingResidues = {}
    fixer.applyMutations([f'ALA-2-{res1}', f'ALA-3-{res2}'], 'A')
    fixer.findMissingAtoms()
    fixer.addMissingAtoms()
    modeller = app.Modeller(fixer.topology, fixer.positions)
    modeller.addHydrogens(forcefield=forcefield, variants=[None, var1, var2, None])
    return modeller.topology, modeller.positions

def createConformations(forcefield, topology, positions, name, charges):
    """Generate the conformations for a molecule and compute the energies and forces."""
    print(f'Generating {name}')
    system = forcefield.createSystem(topology, constraints=None, rigidWater=False)
    integrator = openmm.LangevinMiddleIntegrator(500*unit.kelvin, 1/unit.picosecond, 0.001*unit.picosecond)
    simulation = app.Simulation(topology, system, integrator, openmm.Platform.getPlatformByName('Reference'))
    simulation.context.setPositions(positions)
    simulation.minimizeEnergy()

    # Generate 10 diverse starting points.  Run MD from each one to generate a total
    # of 100 conformations at each temperature.

    mol = convert_to_openff(simulation, charges)
    mol.generate_conformers(n_conformers=10, rms_cutoff=0*ffunit.nanometers)
    assert len(mol.conformers) == 10
    states = []
    for temperature in [100, 300, 500, 1000]:
        for pos in mol.conformers:
            simulation.context.setPositions(pos.m_as(ffunit.nanometers))
            simulation.minimizeEnergy()
            simulation.context.setVelocitiesToTemperature(temperature*unit.kelvin)
            integrator.setTemperature(temperature*unit.kelvin)
            for i in range(10):
                simulation.step(10000)
                state = simulation.context.getState(getPositions=True, getEnergy=True)
                if state.getPotentialEnergy() < 1e4*unit.kilojoules_per_mole:
                    states.append(state)

    # Add the conformations to the Molecule and put the atoms in canonical order.

    mol._conformers = None
    for state in states:
        mol.add_conformer(state.getPositions(asNumpy=True))
    mol = mol.canonical_order_atoms()

    # Compute energies and forces with XTB.

    positions, energies, formation_energies, grads = compute_xtb_for_conformers(mol)
    return positions, energies, formation_energies, grads, mol, name

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

forcefield = app.ForceField('amber14-all.xml')
pdb = app.PDBFile('../dipeptides/disulfide.pdb')
futures = []
with ThreadPoolExecutor(max_workers=os.cpu_count()) as executor:
    futures.append(executor.submit(createConformations, forcefield, pdb.topology, pdb.positions, 'Disulfide', {}))
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
            futures.append(executor.submit(createConformations, forcefield, topology, positions, f'{name1}-{name2}', charges))

# Save the results to the output file.

outputfile = h5py.File('dipeptides.hdf5', 'w')
for future in futures:
    positions, energies, formation_energies, grads, mol, name = future.result()
    save_to_file(outputfile, mol, positions, energies, formation_energies, grads, name)
