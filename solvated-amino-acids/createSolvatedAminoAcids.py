import openmm
import openmm.app as app
import openmm.unit as unit
import numpy as np
import h5py
import os
from collections import namedtuple

def createModel(forcefield, res, var):
    """Create a capped amino acid of the specified type and build a water box around it."""
    import pdbfixer
    fixer = pdbfixer.PDBFixer(filename='ace_ala_nme.pdb')
    fixer.missingResidues = {}
    fixer.applyMutations([f'ALA-2-{res}'], 'A')
    fixer.findMissingAtoms()
    fixer.addMissingAtoms()
    modeller = app.Modeller(fixer.topology, fixer.positions)
    modeller.addHydrogens(forcefield=forcefield, variants=[None, var, None])
    modeller.addSolvent(forcefield, boxSize=(2.2, 2.2, 2.2)*unit.nanometer)
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
    chargedAtoms = []
    for i, atom in enumerate(simulation.topology.atoms()):
        obatom = obmol.GetAtom(i+1)

        # Find the formal charge.  In some cases, the PDB file isn't sufficient to determine
        # which atom the charge goes on.  For example in GLU, nothing in the PDB file distinguishes
        # OE1 and OE2.  In these cases we list both atoms but mark them as not required.  We can identify
        # which atom OpenBabel selected because it will have recorded an implicit H (for negative sites)
        # or a valence of 4 (for positive nitrogens).

        charge = 0
        if atom.index in charges:
            c, required = charges[atom.index]
            if required or (c < 0 and obatom.GetImplicitHCount() == 1) or (c > 0 and obatom.GetTotalValence() > 3):
                charge = c
                chargedAtoms.append(obatom)
        stereo = None
        if facade.HasTetrahedralStereo(obatom.GetId()):
            stereo = ('S', 'R')[facade.GetTetrahedralStereo(obatom.GetId()).GetConfig().winding]
        mol.add_atom(atom.element.atomic_number, charge, obatom.IsAromatic(), stereo, atom.name)
    for i in range(obmol.NumBonds()):
        bond = obmol.GetBond(i)
        order = bond.GetBondOrder()
        if (bond.GetBeginAtom().GetImplicitHCount() != 0 and bond.GetEndAtom() in chargedAtoms) or (bond.GetEndAtom().GetImplicitHCount() != 0 and bond.GetBeginAtom() in chargedAtoms):
            # Workaround for a case that comes up in HIP.  OpenBabel infers an implicit hydrogen
            # on CE1, when instead we want it to infer a +1 charge on ND1, leading to the wrong
            # order for the bond connecting them.
            order += 1
        mol.add_bond(bond.GetBeginAtomIdx()-1, bond.GetEndAtomIdx()-1, order, bond.IsAromatic())

    # Verify that the total charge is correct.

    nb = [f for f in simulation.context.getSystem().getForces() if isinstance(f, openmm.NonbondedForce)][0]
    charge = round(sum(nb.getParticleParameters(i)[0]._value for i in range(nb.getNumParticles())))
    assert charge == mol.total_charge.value_in_unit(unit.elementary_charge)
    return mol

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

    # Simulate for 1 ns, saving a state every 20 ps.

    import mdtraj
    atoms = list(topology.atoms())
    soluteAtomIndex = [atom.index for atom in atoms if atom.residue.chain.index == 0]
    waterAtomIndex = [atom.index for atom in atoms if atom.residue.name == 'HOH' and atom.name == 'O']
    from openmm.app.internal import compiled
    conformations = []
    for i in range(50):
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

    # Create an OpenFF molecule and save the results.

    system = forcefield.createSystem(modeller.topology, constraints=None, rigidWater=False)
    integrator = openmm.LangevinMiddleIntegrator(300*unit.kelvin, 1/unit.picosecond, 0.001*unit.picosecond)
    simulation = app.Simulation(modeller.topology, system, integrator)
    simulation.context.setPositions(conformations[-1])
    mol = convertToOpenFF(simulation, charges)
    saveToFile(outputfile, mol, conformations, name)

def saveToFile(outputfile, mol, conformations, name):
    """Save a molecule and its conformations to a HDF5 file."""
    mol._conformers = None
    for conf in conformations:
        mol.add_conformer(conf)
    mol = mol.canonical_order_atoms()
    smiles = mol.to_smiles(isomeric=True, explicit_hydrogens=True, mapped=True)
    conformations = [c.value_in_unit(unit.nanometers) for c in mol.conformers]
    group = outputfile.create_group(name)
    group.create_dataset('smiles', data=[smiles], dtype=h5py.string_dtype())
    ds = group.create_dataset('conformations', data=np.array(conformations), dtype=np.float32)
    ds.attrs['units'] = 'nanometers'

    # As a sanity check, make sure the SMILES string doesn't have any radicals.

    from rdkit import Chem
    rdmol = Chem.MolFromSmiles(smiles)
    assert all(atom.GetNumRadicalElectrons() == 0 for atom in rdmol.GetAtoms())

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

# Create the molecules.

forcefield = app.ForceField('amber14-all.xml', 'amber14/tip3pfb.xml')
outputfile = h5py.File('solvated-amino-acids.hdf5', 'w')
for res in residues:
    name = res.name if res.variant is None else res.variant
    topology, positions = createModel(forcefield, res.name, res.variant)
    charges = {}
    for atom in topology.atoms():
        if atom.residue.index == 1 and atom.name in res.charges:
            charges[atom.index] = res.charges[atom.name]
    createConformations(outputfile, forcefield, topology, positions, name, charges)
