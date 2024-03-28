import openmm
import openmm.app as app
import openmm.unit as unit
from openff.toolkit.topology import Molecule
from openff.units import unit as ffunit
from rdkit import Chem
from rdkit.Chem import AllChem, rdDetermineBonds
import ase
import ase.md
import ase.md.velocitydistribution
import ase.optimize
from xtb.ase.calculator import XTB
import pdbfixer
import numpy as np
import h5py
import logging
import random
import os

logging.getLogger().setLevel(logging.ERROR)

elements = set(['B', 'Br', 'C', 'Cl', 'F', 'H', 'I', 'N', 'O', 'P', 'S', 'Si'])

def createPeptideSimulation(forcefield, residues):
    """Use PDBFixer to create a peptide with specified residues."""
    fixer = pdbfixer.PDBFixer(filename='../dipeptides/ala_ala.pdb')
    fixer.missingResidues = {(0, 4):residues[3:]+['NME']}
    fixer.applyMutations([f'ALA-2-{residues[0]}', f'ALA-3-{residues[1]}', f'NME-4-{residues[2]}'], 'A')
    fixer.findMissingAtoms()
    fixer.addMissingAtoms()
    fixer.addMissingHydrogens(7.0, forcefield)
    system = forcefield.createSystem(fixer.topology)
    integrator = openmm.LangevinMiddleIntegrator(300*unit.kelvin, 1/unit.picosecond, 0.001*unit.picosecond)
    simulation = app.Simulation(fixer.topology, system, integrator, openmm.Platform.getPlatformByName('Reference'))
    simulation.context.setPositions(fixer.positions)
    return simulation

def createPeptideMolecule(simulation):
    """Create an OpenFF molecule for a peptide."""
    rdmol = Chem.EditableMol(Chem.Mol())
    for atom in simulation.topology.atoms():
        a = Chem.Atom(atom.element.atomic_number)
        a.SetNoImplicit(True)
        rdmol.AddAtom(a)
    for bond in simulation.topology.bonds():
        rdmol.AddBond(bond[0].index, bond[1].index, Chem.BondType.SINGLE)
    rdmol = rdmol.GetMol()
    nonbonded = [f for f in simulation.system.getForces() if isinstance(f, openmm.NonbondedForce)][0]
    charge = round(sum(nonbonded.getParticleParameters(i)[0].value_in_unit(unit.elementary_charge) for i in range(nonbonded.getNumParticles())))
    rdDetermineBonds.DetermineBondOrders(rdmol, charge, embedChiral=False)
    return Molecule(rdmol, allow_undefined_stereo=True)

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

    rdmol = Chem.MolFromSmiles(smiles)
    assert all(atom.GetNumRadicalElectrons() == 0 for atom in rdmol.GetAtoms())

def processLigand(smiles, mol):
    """Generate conformations for a ligand."""
    print(f'Generating {smiles}')
    AllChem.EmbedMultipleConfs(mol, numConfs=10, pruneRmsThresh=0)
    assert mol.GetNumConformers() == 10
    numbers = [a.GetAtomicNum() for a in mol.GetAtoms()]
    charges = [a.GetFormalCharge() for a in mol.GetAtoms()]
    conformations = []
    for conf in mol.GetConformers():
        atoms = ase.Atoms(positions=conf.GetPositions(), numbers=numbers, charges=charges)
        atoms.calc = XTB(method="GFN-FF")
        ase.optimize.LBFGS(atoms, logfile=os.devnull).run(0.001, 20)
        ase.md.velocitydistribution.MaxwellBoltzmannDistribution(atoms, temperature_K=300)
        dyn = ase.md.langevin.Langevin(atoms, 1*ase.units.fs, temperature_K=300, friction=1e-3)
        dyn.run(5000)
        conformations.append(atoms.get_positions()-np.mean(atoms.get_positions(), axis=0))
    return conformations

def processPeptide(forcefield, residues):
    """Generate conformations for a peptide."""
    print(f'Generating {"-".join(residues)}')
    simulation = createPeptideSimulation(forcefield, residues)
    if simulation.system.getNumParticles() > 110:
        raise ValueError('Skipping too large peptide')
    simulation.minimizeEnergy()
    simulation.context.setVelocitiesToTemperature(300*unit.kelvin)
    conformations = []
    for i in range(10):
        simulation.step(50000)
        conformations.append(simulation.context.getState(getPositions=True).getPositions(asNumpy=True).value_in_unit(unit.angstrom))
    return conformations, createPeptideMolecule(simulation)

# Read the list of ligands, filter them to find ones we want to include, and process them.

outputfile = h5py.File('test-dataset.hdf5', 'w')
lines = list(open('../amino-acid-ligand/Components-smiles-oe.smi'))
random.Random(2).shuffle(lines)
numSmall = 0
numLarge = 0
for line in lines:
    smiles, resid, name = line.split('\t')
    resid = resid.upper()
    if resid in app.PDBFile._standardResidues:
        continue
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        continue
    if any(a.GetSymbol() not in elements for a in mol.GetAtoms()):
        continue
    mol = Chem.AddHs(mol)
    if any(a.GetNumRadicalElectrons() != 0 or a.GetIsotope() != 0 for a in mol.GetAtoms()):
        continue
    if mol.GetNumAtoms() >= 40 and mol.GetNumAtoms() <= 50 and numSmall < 200:
        try:
            conformations = processLigand(smiles, mol)
            saveToFile(outputfile, Molecule(mol, allow_undefined_stereo=True), conformations, resid)
            numSmall += 1
        except Exception as ex:
            print(ex)
            pass
    if mol.GetNumAtoms() >= 70 and mol.GetNumAtoms() <= 80 and numLarge < 200:
        try:
            conformations = processLigand(smiles, mol)
            saveToFile(outputfile, Molecule(mol, allow_undefined_stereo=True), conformations, resid)
            numLarge += 1
        except Exception as ex:
            print(ex)
            pass

# Generate random peptides and process them.

forcefield = app.ForceField('amber14/protein.ff14SB.xml')
numPeptides = 0
while numPeptides < 200:
    residues = []
    for j in range(5):
        residues.append(random.choice(pdbfixer.pdbfixer.proteinResidues))
    try:
        conformations, mol = processPeptide(forcefield, residues)
        saveToFile(outputfile, mol, conformations, '-'.join(residues))
        numPeptides += 1
    except Exception as ex:
        print(ex)
        pass
