from qcportal import PortalClient
from collections import defaultdict
from rdkit import Chem
import numpy as np
import h5py
import yaml

# Units for a variety of fields that can be downloaded.

units = {'dft_total_energy': 'hartree',
         'dft_total_gradient': 'hartree/bohr',
         'mbis_charges': 'elementary_charge',
         'mbis_dipoles': 'elementary_charge*bohr',
         'mbis_quadrupoles': 'elementary_charge*bohr^2',
         'mbis_octupoles': 'elementary_charge*bohr^3',
         'scf_dipole': 'elementary_charge*bohr',
         'scf_quadrupole': 'elementary_charge*bohr^2'}

# Reference energies computed with Psi4 1.5 wB97M-D3BJ/def2-TZVPPD.

atom_energy = {'B': {-1: -24.677421752684776, 0: -24.671520535482145, 1: -24.364648707125294},
               'Br': {-1: -2574.2451510945853, 0: -2574.1167240829964},
               'C': {-1: -37.91424135791358, 0: -37.87264507233593, 1: -37.45349214963933},
               'Ca': {2: -676.9528465198214},
               'Cl': {-1: -460.3350243496703, 0: -460.1988762285739},
               'F': {-1: -99.91298732343974, 0: -99.78611622985483},
               'H': {-1: -0.5027370838721259, 0: -0.4987605100487531, 1: 0.0},
               'I': {-1: -297.8813829975981, 0: -297.76228914445625},
               'K': {1: -599.8025677513111},
               'Li': {1: -7.285254714046546},
               'Mg': {2: -199.2688420040449},
               'N': {-1: -54.602291095426494, 0: -54.62327513368922, 1: -54.08594142587869},
               'Na': {1: -162.11366478783253},
               'O': {-1: -75.17101657391741, 0: -75.11317840410095, 1: -74.60241514396725},
               'P': {0: -341.3059197024934, 1: -340.9258392474849},
               'S': {-1: -398.2405387031612, 0: -398.1599636677874, 1: -397.7746615977658},
               'Si': {-1: -289.4540686037408, 0: -289.4131352299586, 1: -289.1189404777897}}

default_charge = {}
for symbol in atom_energy:
    energies = [(energy, charge) for charge, energy in atom_energy[symbol].items()]
    default_charge[symbol] = sorted(energies)[0][1]

# Given a SMILES string, compute the reference energy of the atoms when fully separated.

def compute_reference_energy(smiles):
    rdmol = Chem.MolFromSmiles(smiles, sanitize=False)
    total_charge = sum(atom.GetFormalCharge() for atom in rdmol.GetAtoms())
    symbol = [atom.GetSymbol() for atom in rdmol.GetAtoms()]
    charge = [default_charge[s] for s in symbol]
    delta = np.sign(total_charge-sum(charge))
    while delta != 0:
        best_index = -1
        best_energy = None
        for i in range(len(symbol)):
            s = symbol[i]
            e = atom_energy[s]
            new_charge = charge[i]+delta
            if new_charge in e:
                if best_index == -1 or e[new_charge]-e[charge[i]] < best_energy:
                    best_index = i
                    best_energy = e[new_charge]-e[charge[i]]
        charge[best_index] += delta
        delta = np.sign(total_charge - sum(charge))
    return sum(atom_energy[s][c] for s, c in zip(symbol, charge))

# Assemble the value of a data field from the record properties.

def get_data_value(name, qcvars):
    altname = name.replace(' ', '_')
    data = []
    for vars in qcvars:
        if name in vars:
            value = vars[name]
        else:
            value = vars[altname]
        if 'gradient' in name:
            value = np.reshape(value, [-1, 3])
        elif 'charges' in name:
            value = np.reshape(value, [-1, 1])
        elif 'dipoles' in name:
            value = np.reshape(value, [-1, 3])
        elif 'quadrupoles' in name:
            value = np.reshape(value, [-1, 3, 3])
        elif 'octupoles' in name:
            value = np.reshape(value, [-1, 3, 3, 3])
        elif altname == 'scf_quadrupole':
            value = np.reshape(value, [3, 3])
        elif 'indices' in name:
            atoms = int(np.sqrt(len(value)))
            value = np.reshape(value, [atoms, atoms])
        data.append(value)
    return data

# Process the configuration file and download data.

with open('config.yaml') as input:
    config = yaml.safe_load(input.read())
if 'max_force' in config:
    max_force = float(config['max_force'])
else:
    max_force = None
client = PortalClient('https://ml.qcarchive.molssi.org')
outputfile = h5py.File('SPICE.hdf5', 'w')
for subset in config['subsets']:
    # Download the next subset.

    print('Processing', subset)
    dataset = client.get_dataset('singlepoint', subset)
    specifications = [name for name in dataset.specifications if dataset.specifications[name].specification.method == 'wb97m-d3bj']
    while True:
        try:
            recs = list(dataset.iterate_records(specification_names=specifications))
            break
        except Exception as error:
            print(error)
            print('Retrying')
    recs_by_name = defaultdict(list)
    for e, s, r in recs:
        if r is not None and r.status == 'complete':
            name = e[:e.rfind('-')]
            recs_by_name[name].append(r)
    all_molecules = client.get_molecules([r.molecule_id for e, s, r in recs])
    mols_by_id = dict((m.id, m) for m in all_molecules)

    # Add the data to the HDF5 file.

    for name in recs_by_name:
        group_recs = recs_by_name[name]
        molecules = [mols_by_id[r.molecule_id] for r in group_recs]
        qcvars = [r.properties for r in group_recs]
        smiles = molecules[0].extras['canonical_isomeric_explicit_hydrogen_mapped_smiles']
        ref_energy = compute_reference_energy(smiles)
        name = name.replace('/', '')  # Remove stereochemistry markers that h5py interprets as path separators
        group = outputfile.create_group(name)
        group.create_dataset('subset', data=[subset], dtype=h5py.string_dtype())
        group.create_dataset('smiles', data=[smiles], dtype=h5py.string_dtype())
        group.create_dataset("atomic_numbers", data=molecules[0].atomic_numbers, dtype=np.int16)
        if max_force is not None:
            force = np.array([vars['dft total gradient'] for vars in qcvars])
            samples = [i for i in range(len(molecules)) if np.max(np.abs(force[i])) <= max_force]
            molecules = [molecules[i] for i in samples]
            qcvars = [qcvars[i] for i in samples]
        ds = group.create_dataset('conformations', data=np.array([m.geometry for m in molecules]), dtype=np.float32)
        ds.attrs['units'] = 'bohr'
        ds = group.create_dataset('formation_energy', data=np.array([vars['dft total energy']-ref_energy for vars in qcvars]), dtype=np.float64)
        ds.attrs['units'] = 'hartree'
        for value in config['values']:
            altvalue = value.replace(' ', '_')
            key = altvalue.lower()
            try:
                data = get_data_value(value, qcvars)
                dtype = np.float64 if 'energy' in key else np.float32
                ds = group.create_dataset(key, data=np.array(data), dtype=dtype)
                if key in units:
                    ds.attrs['units'] = units[key]
            except KeyError:
                pass
