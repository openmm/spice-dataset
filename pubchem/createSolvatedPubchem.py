import openmm
import openmm.app as app
import openmm.unit as unit
from openff.toolkit.typing.engines.smirnoff import ForceField
from openff.toolkit.topology import Molecule, Topology
from openff.units import unit as ffunit
from openmmforcefields.generators import SMIRNOFFTemplateGenerator
import numpy as np
import h5py
import sys

def topologyToSmiles(topology):
    import copy
    smiles = []
    for top_mol in topology.topology_molecules:
        ref_mol_copy = copy.deepcopy(top_mol.reference_molecule)
        ref_mol_copy.properties["atom_map"] = dict()
        for atom in top_mol.atoms:
            ref_mol_copy.properties["atom_map"][
                atom.atom.molecule_atom_index
            ] = atom.topology_atom_index
        smiles.append(ref_mol_copy.to_smiles(isomeric=True, explicit_hydrogens=True, mapped=True))
    return ".".join(smiles)

def createConformations(outputfile, smiles, sid):
    """Generate the conformations for a molecule and save them to disk."""
    print(f'Generating {index}: {smiles}')
    mol = Molecule.from_smiles(smiles, allow_undefined_stereo=True)
    mol.generate_conformers(n_conformers=1)
    fftop = Topology()
    fftop.add_molecule(mol)
    mmtop = fftop.to_openmm()
    smirnoff = SMIRNOFFTemplateGenerator(molecules=mol)
    forcefield = app.ForceField('amber14/tip3pfb.xml')
    forcefield.registerTemplateGenerator(smirnoff.generator)
    modeller = app.Modeller(mmtop, mol.conformers[0].m_as(ffunit.nanometer))
    try:
        modeller.addSolvent(forcefield, boxSize=(2.2, 2.2, 2.2)*unit.nanometer)
        mmtop = modeller.topology
        system = forcefield.createSystem(mmtop, nonbondedMethod=app.CutoffPeriodic, constraints=None, rigidWater=True)
    except:
        print('  failed to parametrize')
        return
    system.addForce(openmm.MonteCarloBarostat(1*unit.bar, 300*unit.kelvin))
    integrator = openmm.LangevinMiddleIntegrator(300*unit.kelvin, 1/unit.picosecond, 0.001*unit.picosecond)
    simulation = app.Simulation(modeller.topology, system, integrator)
    simulation.context.setPositions(modeller.positions)

    # Equilibrate for 100 ps.

    simulation.minimizeEnergy()
    simulation.context.setVelocitiesToTemperature(300*unit.kelvin)
    simulation.step(100000)

    # Simulate for 200 ps, saving a state every 20 ps.

    import mdtraj
    atoms = list(mmtop.atoms())
    soluteAtomIndex = [atom.index for atom in atoms if atom.residue.chain.index == 0]
    waterAtomIndex = [atom.index for atom in atoms if atom.residue.name == 'HOH' and atom.name == 'O']
    from openmm.app.internal import compiled
    conformations = []
    for i in range(10):
        simulation.step(20000)
        state = simulation.context.getState(getPositions=True, getEnergy=True)
        assert state.getPotentialEnergy() < 0*unit.kilojoules_per_mole

        # Identify the 20 water molecules closest to the solute.

        periodicDistance = compiled.periodicDistance(state.getPeriodicBoxVectors().value_in_unit(unit.nanometer))
        pos = state.getPositions().value_in_unit(unit.nanometer)
        distances = [min(periodicDistance(pos[i], pos[j]) for j in soluteAtomIndex) for i in waterAtomIndex]
        sortedIndices = np.argsort(distances)
        waterToKeep = set(atoms[waterAtomIndex[i]].residue.index for i in sortedIndices[:20])

        # Remove everything else from the Topology and coordinates.

        modeller = app.Modeller(mmtop, state.getPositions())
        toDelete = [res for res in mmtop.residues() if res.chain.index > 0 and res.index not in waterToKeep]
        modeller.delete(toDelete)

        # Center the solute with the water around it.

        mdtop = mdtraj.Topology.from_openmm(modeller.topology)
        xyz = np.array([modeller.positions.value_in_unit(unit.nanometer)])
        traj = mdtraj.Trajectory(xyz, mdtop)
        traj.unitcell_vectors = np.array([state.getPeriodicBoxVectors().value_in_unit(unit.nanometer)])
        mdmols = mdtop.find_molecules()
        traj = traj.image_molecules(anchor_molecules=[mdmols[0]], other_molecules=mdmols[1:], make_whole=False)
        conformations.append(traj.xyz[0])

    # Order the atoms within the PubChem molecule.

    mol._conformers = None
    for c in conformations:
        mol.add_conformer(c[:mol.n_atoms]*ffunit.nanometers)
    mol = mol.canonical_order_atoms()

    # Generate the SMILES string.

    smiles = mol.to_smiles(isomeric=True, explicit_hydrogens=True, mapped=True)
    for i in range(20):
        start = mol.n_atoms+3*i+1
        smiles += f'.[O:{start}]([H:{start+1}])[H:{start+2}]'

    # Generate the complete conformations and save them to the file.

    for i in range(len(conformations)):
        conformations[i][:mol.n_atoms] = mol.conformers[i].m_as(ffunit.nanometers)
    conformations = [c-np.average(c, axis=0) for c in conformations]
    group = outputfile.create_group(f'solvated {sid}')
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
outputfile = h5py.File(f'solvated-pubchem-{first}-{last}.hdf5', 'w')
index = 0
for line in open('sorted.txt'):
    index += 1
    if index < first:
        continue
    if index > last:
        break
    sid, smiles = line.split()
    try:
        createConformations(outputfile, smiles, sid)
    except:
        print('  failed to generate conformations')
