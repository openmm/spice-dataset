import openmm
import openmm.app as app
import openmm.unit as unit
from openff.toolkit.topology import Molecule
from openff.units import unit as ffunit
from openmmforcefields.generators import SMIRNOFFTemplateGenerator
from rdkit import Chem
import pdbfixer
import numpy as np
import h5py
from collections import defaultdict
from io import StringIO
import logging
import os
import tempfile
import urllib.request

logging.getLogger().setLevel(logging.ERROR)

elements = set(['B', 'Br', 'C', 'Cl', 'F', 'H', 'I', 'N', 'O', 'P', 'S', 'Si'])
aminoAcidSmiles = {
    'ALA': '[H][N]([C](=[O])[C@]([H])([N]([H])[C](=[O])[C]([H])([H])[H])[C]([H])([H])[H])[C]([H])([H])[H]',
    'ARG': '[H][N]([H])[C]([N]([H])[C]([H])([H])[C]([H])([H])[C]([H])([H])[C@]([H])([C](=[O])[N]([H])[C]([H])([H])[H])[N]([H])[C](=[O])[C]([H])([H])[H])=[N+]([H])[H]',
    'ASN': '[H][N]([H])[C](=[O])[C]([H])([H])[C@]([H])([C](=[O])[N]([H])[C]([H])([H])[H])[N]([H])[C](=[O])[C]([H])([H])[H]',
    'ASP': '[H][N]([C](=[O])[C@]([H])([N]([H])[C](=[O])[C]([H])([H])[H])[C]([H])([H])[C](=[O])[O-])[C]([H])([H])[H]',
    'CYS': '[H][S][C]([H])([H])[C@@]([H])([C](=[O])[N]([H])[C]([H])([H])[H])[N]([H])[C](=[O])[C]([H])([H])[H]',
    'GLN': '[H][N]([H])[C](=[O])[C]([H])([H])[C]([H])([H])[C@]([H])([C](=[O])[N]([H])[C]([H])([H])[H])[N]([H])[C](=[O])[C]([H])([H])[H]',
    'GLU': '[H][N]([C](=[O])[C@]([H])([N]([H])[C](=[O])[C]([H])([H])[H])[C]([H])([H])[C]([H])([H])[C](=[O])[O-])[C]([H])([H])[H]',
    'GLY': '[H][N]([C](=[O])[C]([H])([H])[N]([H])[C](=[O])[C]([H])([H])[H])[C]([H])([H])[H]',
    'HIS': '[H][C]1=[C]([C]([H])([H])[C@]([H])([C](=[O])[N]([H])[C]([H])([H])[H])[N]([H])[C](=[O])[C]([H])([H])[H])[N+]([H])=[C]([H])[N]1[H]',
    'ILE': '[H][N]([C](=[O])[C@]([H])([N]([H])[C](=[O])[C]([H])([H])[H])[C@]([H])([C]([H])([H])[H])[C]([H])([H])[C]([H])([H])[H])[C]([H])([H])[H]',
    'LEU': '[H][N]([C](=[O])[C@]([H])([N]([H])[C](=[O])[C]([H])([H])[H])[C]([H])([H])[C]([H])([C]([H])([H])[H])[C]([H])([H])[H])[C]([H])([H])[H]',
    'LYS': '[H][N]([C](=[O])[C@]([H])([N]([H])[C](=[O])[C]([H])([H])[H])[C]([H])([H])[C]([H])([H])[C]([H])([H])[C]([H])([H])[N+]([H])([H])[H])[C]([H])([H])[H]',
    'MET': '[H][N]([C](=[O])[C@]([H])([N]([H])[C](=[O])[C]([H])([H])[H])[C]([H])([H])[C]([H])([H])[S][C]([H])([H])[H])[C]([H])([H])[H]',
    'PHE': '[H][c]1[c]([H])[c]([H])[c]([C]([H])([H])[C@]([H])([C](=[O])[N]([H])[C]([H])([H])[H])[N]([H])[C](=[O])[C]([H])([H])[H])[c]([H])[c]1[H]',
    'PRO': '[H][N]([C](=[O])[C@]1([H])[N]([C](=[O])[C]([H])([H])[H])[C]([H])([H])[C]([H])([H])[C]1([H])[H])[C]([H])([H])[H]',
    'SER': '[H][O][C]([H])([H])[C@]([H])([C](=[O])[N]([H])[C]([H])([H])[H])[N]([H])[C](=[O])[C]([H])([H])[H]',
    'THR': '[H][O][C@]([H])([C]([H])([H])[H])[C@]([H])([C](=[O])[N]([H])[C]([H])([H])[H])[N]([H])[C](=[O])[C]([H])([H])[H]',
    'TRP': '[H][C]1=[C]([C]([H])([H])[C@]([H])([C](=[O])[N]([H])[C]([H])([H])[H])[N]([H])[C](=[O])[C]([H])([H])[H])[c]2[c]([H])[c]([H])[c]([H])[c]([H])[c]2[N]1[H]',
    'TYR': '[H][O][c]1[c]([H])[c]([H])[c]([C]([H])([H])[C@]([H])([C](=[O])[N]([H])[C]([H])([H])[H])[N]([H])[C](=[O])[C]([H])([H])[H])[c]([H])[c]1[H]',
    'VAL': '[H][N]([C](=[O])[C@]([H])([N]([H])[C](=[O])[C]([H])([H])[H])[C]([H])([C]([H])([H])[H])[C]([H])([H])[H])[C]([H])([H])[H]',
}

def createOpenFFMolecule(topology, aaName, ligandPdbPath, smiles):
    """Create an OpenForceField Molecule object combining an amino acid with a ligand."""

    # First create a Molecule for the amino acid.  We can create it from a SMILES string.  This gets
    # all the chemistry right (bond orders, partial charges, etc.), but the atoms are in an unknown
    # order.  So we build a second Molecule from the Topology and have OpenFF work out the mapping
    # between the two.

    mol = Molecule.from_smiles(aminoAcidSmiles[aaName], allow_undefined_stereo=True)
    mol2 = Molecule()
    for atom in topology.atoms():
        mol2.add_atom(atom.element.atomic_number, 0, False, name=atom.name)
    for a1, a2 in topology.bonds():
        mol2.add_bond(a1.index, a2.index, 1, False)
    isomorphic, mapping = Molecule.are_isomorphic(mol, mol2,
        return_atom_map=True,
        aromatic_matching=False,
        formal_charge_matching=False,
        bond_order_matching=False,
        atom_stereochemistry_matching=False,
        bond_stereochemistry_matching=False,
    )
    mol = mol.remap(mapping)

    # Create another Molecule for the ligand and merge them together.

    ligand = Molecule.from_pdb_and_smiles(ligandPdbPath, smiles, allow_undefined_stereo=True)
    offset = topology.getNumAtoms()
    for atom in ligand.atoms:
        mol.add_atom(atom.atomic_number, atom.formal_charge, atom.is_aromatic, atom.stereochemistry, atom.name)
    for bond in ligand.bonds:
        mol.add_bond(bond.atom1_index+offset, bond.atom2_index+offset, bond.bond_order, bond.is_aromatic)
    return mol

def createLigandModel(pdb, ligand, ligandPdb):
    """Create an OpenMM model for just the ligand."""
    # Delete everything except the ligand.

    modeller = app.Modeller(pdb.topology, pdb.positions)
    modeller.delete([r for r in pdb.topology.residues() if r != ligand])

    # The PDB file may be missing hydrogens.  Modeller can add them, but since it's a nonstandard
    # reside, we first need to provide a definition for it.  We can get the hydrogens, atom names,
    # and bonds from the ligand PDB file provided in the Chemical Components Dictionary.

    data = app.Modeller._ResidueData(ligand.id)
    app.Modeller._residueHydrogens[ligand.name] = data
    for a1, a2, in ligandPdb.topology.bonds():
        if a1.element == app.element.hydrogen:
            a1, a2 = a2, a1
        if a2.element == app.element.hydrogen:
            data.hydrogens.append(app.Modeller._Hydrogen(a2.name, a1.name, np.inf, None, None))
    modeller.addHydrogens()
    return modeller

def findNeighbors(test_atoms, test_positions, residue):
    """Find all amino acids within a cutoff distance of the ligand."""
    mindist2 = np.ones(test_positions.shape[0])*np.inf
    for atom in residue.atoms():
        delta = test_positions-positions[atom.index]
        mindist2 = np.minimum(mindist2, np.sum(delta*delta, axis=1))
    mindist = np.sqrt(mindist2)
    neighbors = set()
    for i in range(len(test_atoms)):
        if mindist[i] < 0.4:
            neighbors.add(test_atoms[i].residue)
    return neighbors

def createModel(pdb, ligand, neighbor, ligandModel, forcefield, ligandPdbPath, smiles):
    """Create an OpenMM model for the ligand and amino acid."""
    # Delete everything except the ligand and amino acid.

    modeller = app.Modeller(pdb.topology, pdb.positions)
    modeller.delete([r for r in pdb.topology.residues() if r not in (ligand, neighbor)])

    # Use PDBFixer to cap the amino acid.  Having the ligand present at this point is
    # important, so it can prevent the caps from overlapping with it.

    output = StringIO()
    app.PDBFile.writeFile(modeller.topology, modeller.positions, output)
    fixer = pdbfixer.PDBFixer(pdbfile=StringIO(output.getvalue()))
    aa = [r for r in fixer.topology.residues() if r.name == neighbor.name][0]
    fixer.missingResidues = {(aa.chain.index, 0): ['ACE'], (aa.chain.index, 1): ['NME']}
    fixer.findMissingAtoms()
    fixer.addMissingAtoms()

    # Delete the ligand and add hydrogens to the amino acid.

    modeller = app.Modeller(fixer.topology, fixer.positions)
    modeller.delete([r for r in modeller.topology.residues() if r.name == ligand.name])
    variants = ['HIP' if r.name == 'HIS' else None for r in modeller.topology.residues()]
    modeller.addHydrogens(forcefield, variants=variants)

    # Create both a Molecule and a Topology for the combined system (ligand and amino acid).

    mol = createOpenFFMolecule(modeller.topology, neighbor.name, ligandPdbPath, smiles)
    modeller.add(ligandModel.topology, ligandModel.positions)
    return modeller, mol

def createForceField(smiles):
    """Create a ForceField that can be used to parametrize the amino acid and ligand."""
    molecule = Molecule.from_smiles(smiles, allow_undefined_stereo=True)
    smirnoff = SMIRNOFFTemplateGenerator(molecules=molecule)
    forcefield = app.ForceField('amber14/protein.ff14SB.xml')
    forcefield.registerTemplateGenerator(smirnoff.generator)
    return forcefield

def createSystem(model, forcefield):
    """Create an OpenMM System, including harmonic restraints on heavy atoms."""
    res = list(model.topology.residues())[3]
    system = forcefield.createSystem(model.topology)
    restraint = openmm.CustomExternalForce('1000*((x-x0)^2+(y-y0)^2+(z-z0)^2)')
    restraint.addPerParticleParameter('x0')
    restraint.addPerParticleParameter('y0')
    restraint.addPerParticleParameter('z0')
    for atom in model.topology.atoms():
        if atom.element != app.element.hydrogen and atom.residue.name not in ('ACE', 'NME'):
            restraint.addParticle(atom.index, model.positions[atom.index])
    system.addForce(restraint)
    return system

def createConformation(model, system):
    """Create a single conformation for the model."""
    integrator = openmm.LangevinMiddleIntegrator(300*unit.kelvin, 10/unit.picosecond, 0.001*unit.picosecond)
    simulation = app.Simulation(model.topology, system, integrator)
    simulation.context.setPositions(model.positions)
    simulation.minimizeEnergy()
    return simulation.context.getState(getPositions=True).getPositions(asNumpy=True).value_in_unit(unit.nanometer)

def saveConformations(outputfile, conformationMap, molMap):
    """Save all conformations involving a ligand to the output file."""
    for key in conformationMap:
        mol = molMap[key]
        ligandName, aaName = key
        mol._conformers = None
        for conf in conformationMap[key]:
            mol.add_conformer(conf*ffunit.nanometer)
        mol = mol.canonical_order_atoms()
        try:
            smiles = mol.to_smiles(isomeric=True, explicit_hydrogens=True, mapped=True)
        except:
            print('  exception generating canonical SMILES')
            return
        conformations = [c.m_as(ffunit.nanometers) for c in mol.conformers]
        conformations = [c-np.average(c, axis=0) for c in conformations]
        group = outputfile.create_group(f'{aaName} {ligandName}')
        group.create_dataset('smiles', data=[smiles], dtype=h5py.string_dtype())
        ds = group.create_dataset('conformations', data=np.array(conformations), dtype=np.float32)
        ds.attrs['units'] = 'nanometers'

        # As a sanity check, make sure the SMILES string doesn't have any radicals.

        from rdkit import Chem
        rdmol = Chem.MolFromSmiles(smiles)
        assert all(atom.GetNumRadicalElectrons() == 0 for atom in rdmol.GetAtoms())

ligandPDB = {}
for line in open('cc-to-pdb.tdd'):
    ligandPDB[line[:3].upper()] = line[4:8].upper()

outputfile = h5py.File(f'amino-acid-ligand.hdf5', 'w')

for line in open('Components-smiles-oe.smi'):
    smiles, resid, name = line.split('\t')
    resid = resid.upper()
    if resid in app.PDBFile._standardResidues:
        continue
    if resid not in ligandPDB:
        continue
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        continue
    if any(a.GetSymbol() not in elements for a in mol.GetAtoms()):
        continue
    mol = Chem.AddHs(mol)
    if any(a.GetNumRadicalElectrons() != 0 or a.GetIsotope() != 0 for a in mol.GetAtoms()):
        continue
    if mol.GetNumAtoms() < 5 or mol.GetNumAtoms() > 50:
        continue
    print(f'Generating {smiles}')
    with tempfile.TemporaryDirectory() as tempdir:
        pdbPath = f'{tempdir}/model.pdb'
        ligandPdbPath = f'{tempdir}/ligand.pdb'
        urllib.request.urlretrieve(f'https://files.rcsb.org/download/{ligandPDB[resid]}.pdb', pdbPath)
        urllib.request.urlretrieve(f'http://ligand-expo.rcsb.org/reports/{resid[0]}/{resid}/{resid}_model.pdb', ligandPdbPath)
        pdb = app.PDBFile(pdbPath)
        ligandPdb = app.PDBFile(ligandPdbPath)
        test_atoms = [a for a in pdb.topology.atoms() if a.residue.name in aminoAcidSmiles.keys()]
        positions = pdb.getPositions(asNumpy=True).value_in_unit(unit.nanometer)
        test_positions = positions[[a.index for a in test_atoms]]
        forcefield = createForceField(smiles)
        conformations = defaultdict(list)
        mols = {}
        for residue in pdb.topology.residues():
            if residue.name == resid:
                ligandModel = createLigandModel(pdb, residue, ligandPdb)
                neighbors = findNeighbors(test_atoms, test_positions, residue)
                for neighbor in neighbors:
                    try:
                        model, mol = createModel(pdb, residue, neighbor, ligandModel, forcefield, ligandPdbPath, smiles)
                        system = createSystem(model, forcefield)
                        conformations[(residue.name, neighbor.name)].append(createConformation(model, system))
                        mols[(residue.name, neighbor.name)] = mol
                    except:
                        # This happens when model.pdb and ligand.pdb don't have identical sets of heavy atoms.
                        pass
    if len(mols) > 0:
        saveConformations(outputfile, conformations, mols)