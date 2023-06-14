import openmm
import openmm.app as app
import openmm.unit as unit
import h5py
import os
from utils import compute_xtb_for_conformers, save_to_file
from rdkit import Chem
from openff.units import unit as ffunit
from openff.toolkit.typing.engines.smirnoff import ForceField
from openff.toolkit.topology import Molecule, Topology
from concurrent.futures import ThreadPoolExecutor

def createConformations(topology, system, mol, name):
    """Generate the conformations for a molecule and compute the energies and forces."""
    print(f'Generating {name}')
    try:
        system.addForce(openmm.CMMotionRemover())
        integrator = openmm.LangevinMiddleIntegrator(500*unit.kelvin, 1/unit.picosecond, 0.001*unit.picosecond)
        simulation = app.Simulation(topology, system, integrator, openmm.Platform.getPlatformByName('Reference'))

        # Generate 10 diverse starting points.  Run MD from each one to generate a total
        # of 100 conformations at each temperature.

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
    except Exception as ex:
        print(name, ex)
        return None, None, None, None, mol, name

# Create the conformations for the molecules.

forcefield = ForceField('openff_unconstrained-2.0.0.offxml')
outputfile = h5py.File('des-monomers.hdf5', 'w')
futures = []
with ThreadPoolExecutor(max_workers=os.cpu_count()) as executor:
    for filename in os.listdir('../des370k/SDFS'):
        if not filename.endswith('.sdf'):
            continue
        smiles = filename[:-4]
        supp = Chem.SDMolSupplier(f'../des370k/SDFS/{filename}', sanitize=False, removeHs=False)
        rdmol = list(supp)[0]
        if rdmol.GetNumAtoms() > 1:
            try:
                mol = Molecule.from_rdkit(rdmol, allow_undefined_stereo=True, hydrogens_are_explicit=True)
                topology = Topology.from_molecules([mol])
                mmTopology = topology.to_openmm()
                system = forcefield.create_openmm_system(topology)
                futures.append(executor.submit(createConformations, mmTopology, system, mol, smiles))
            except:
                print('  failed to parametrize')

# Save the results to the output file.

for future in futures:
    positions, energies, formation_energies, grads, mol, name = future.result()
    if positions != None:
        save_to_file(outputfile, mol, positions, energies, formation_energies, grads, name, 'DES Monomers')
