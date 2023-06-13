import numpy as np
import openmm
import h5py
from openmm.app import element, PDBFile
from openff.units import unit as ffunit
from xtb.interface import Calculator, Param, XTBException
from xtb.libxtb import VERBOSITY_MUTED
from collections import defaultdict

typeDict = {('Br', -1): 0, ('Br', 0): 1, ('C', -1): 2, ('C', 0): 3, ('C', 1): 4, ('Ca', 2): 5, ('Cl', -1): 6,
            ('Cl', 0): 7, ('F', -1): 8, ('F', 0): 9, ('H', 0): 10, ('I', -1): 11, ('I', 0): 12, ('K', 1): 13,
            ('Li', 1): 14, ('Mg', 2): 15, ('N', -1): 16, ('N', 0): 17, ('N', 1): 18, ('Na', 1): 19, ('O', -1): 20,
            ('O', 0): 21, ('O', 1): 22, ('P', 0): 23, ('P', 1): 24, ('S', -1): 25, ('S', 0): 26, ('S', 1): 27}

atom_energy = defaultdict(dict)
for symbol, charge in typeDict:
    for m in [0, 1, 2, 3]:
        try:
            calc = Calculator(Param.GFN2xTB, np.array([element.get_by_symbol(symbol).atomic_number]), np.array([[0, 0, 0]]), charge, m)
            calc.set_verbosity(VERBOSITY_MUTED)
            res = calc.singlepoint()
            atom_energy[symbol][charge] = res.get_energy()
            break
        except:
            pass

default_charge = {}
for symbol in atom_energy:
    energies = [(energy, charge) for charge, energy in atom_energy[symbol].items()]
    default_charge[symbol] = sorted(energies)[0][1]

def compute_reference_energy(symbols, total_charge):
    """Compute the reference energy of a molecule's atoms when fully separated."""
    charge = [default_charge[s] for s in symbols]
    delta = np.sign(total_charge-sum(charge))
    while delta != 0:
        best_index = -1
        best_energy = None
        for i in range(len(symbols)):
            s = symbols[i]
            e = atom_energy[s]
            new_charge = charge[i]+delta
            if new_charge in e:
                if best_index == -1 or e[new_charge]-e[charge[i]] < best_energy:
                    best_index = i
                    best_energy = e[new_charge]-e[charge[i]]
        charge[best_index] += delta
        delta = np.sign(total_charge - sum(charge))
    return sum([atom_energy[s][c] for s, c in zip(symbols, charge)])

def compute_xtb(symbols, positions, total_charge):
    """Compute the energy and forces for a molecule with GFN2-xTB."""
    numbers = np.array([element.get_by_symbol(s).atomic_number for s in symbols])
    calc = Calculator(Param.GFN2xTB, numbers, positions, total_charge)
    calc.set_verbosity(VERBOSITY_MUTED)
    res = calc.singlepoint()
    energy = res.get_energy()
    grad = res.get_gradient()
    formation_energy = energy-compute_reference_energy(symbols, total_charge)
    return energy, formation_energy, grad

def compute_xtb_for_conformers(mol):
    """Compute the energy and forces for every conformer of a molecule with GFN2-xTB."""
    symbols = np.array([atom.symbol for atom in mol.atoms])
    positions = []
    energies = []
    formation_energies = []
    grads = []
    for conf in mol.conformers:
        pos = conf.m_as(ffunit.bohr)
        try:
            energy, formation_energy, grad = compute_xtb(symbols, pos, mol.total_charge.m)
            positions.append(pos)
            energies.append(energy)
            formation_energies.append(formation_energy)
            grads.append(grad)
        except XTBException as ex:
            print(ex)
    return positions, energies, formation_energies, grads

def convert_to_openff(simulation, charges):
    """Create an OpenForceField Molecule object from an OpenMM simulation."""
    # Use Open Babel to identify bond orders and stereochemistry.

    from openbabel import openbabel
    obConversion = openbabel.OBConversion()
    obConversion.SetInFormat("pdb")
    from io import StringIO
    io = StringIO()
    positions = simulation.context.getState(getPositions=True).getPositions()
    PDBFile.writeFile(simulation.topology, positions, io)
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
    assert charge == mol.total_charge.m
    return mol

def save_to_file(outputfile, mol, positions, energies, formation_energies, grads, name):
    """Save a molecule and its conformations to a HDF5 file."""
    print(f'Saving {name}')
    smiles = mol.to_smiles(isomeric=True, explicit_hydrogens=True, mapped=True)
    group = outputfile.create_group(name)
    group.create_dataset('subset', data=['Dipeptides'], dtype=h5py.string_dtype())
    group.create_dataset('smiles', data=[smiles], dtype=h5py.string_dtype())
    group.create_dataset('atomic_numbers', data=[atom.atomic_number for atom in mol.atoms], dtype=np.int16)
    ds = group.create_dataset('conformations', data=np.array(positions), dtype=np.float32)
    ds.attrs['units'] = 'bohr'
    ds = group.create_dataset('formation_energy', data=np.array(formation_energies), dtype=np.float32)
    ds.attrs['units'] = 'hartree'
    ds = group.create_dataset('xtb_total_energy', data=np.array(energies), dtype=np.float32)
    ds.attrs['units'] = 'hartree'
    ds = group.create_dataset('xtb_total_gradient', data=np.array(grads), dtype=np.float32)
    ds.attrs['units'] = 'hartree/bohr'

    # As a sanity check, make sure the SMILES string doesn't have any radicals.

    from rdkit import Chem
    rdmol = Chem.MolFromSmiles(smiles)
    assert all(atom.GetNumRadicalElectrons() == 0 for atom in rdmol.GetAtoms())
