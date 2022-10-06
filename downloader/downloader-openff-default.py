from qcportal import FractalClient
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

atom_energy = {'Br': {-1: -2574.2451510945853, 0: -2574.1167240829964},
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
               'S': {-1: -398.2405387031612, 0: -398.1599636677874, 1: -397.7746615977658}}
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

# Process the configuration file and download data.

with open('config.yaml') as input:
    config = yaml.safe_load(input.read())
if 'max_force' in config:
    max_force = float(config['max_force'])
else:
    max_force = None
client = FractalClient()
outputfile = h5py.File('SPICE.hdf5', 'w')
for subset in config['subsets']:
    # Download the next subset.

    print('Processing', subset)
    ds = client.get_collection('Dataset', subset)
    all_molecules = ds.get_molecules()
    for row in ds.list_records().iloc:
        spec = row.to_dict()
        if spec['method'] == 'b3lyp-d3bj':
            recs = ds.get_records(method=spec['method'], basis=spec['basis'], program=spec['program'], keywords=spec['keywords'])
            break
    recs_by_name = defaultdict(list)
    mols_by_name = defaultdict(list)
    for i in range(len(recs[0])):
        rec_d3bj  = recs[0].iloc[i].record
        rec_b3lyp = recs[1].iloc[i].record
        if rec_d3bj is not None and rec_b3lyp is not None and rec_d3bj.status == 'COMPLETE' and rec_b3lyp.status == 'COMPLETE':
            assert recs[0].index[i] == recs[1].index[i], print("#{}: index does not match".format(i))
            index = recs[0].index[i]
            name = index[:index.rfind('-')]
            #print(recs[1].iloc[i].name, "\t", name)
            recs_by_name[name].append([rec_d3bj, rec_b3lyp])
            mols_by_name[name].append(all_molecules.loc[index][0])

    # Add the data to the HDF5 file.

    for name in recs_by_name:
        group_recs = recs_by_name[name]
        molecules = mols_by_name[name]
        qcvars_d3bj  = [r[0].extras['qcvars'] for r in group_recs]
        qcvars_b3lyp = [r[1].extras['qcvars'] for r in group_recs]
        smiles = molecules[0].extras['canonical_isomeric_explicit_hydrogen_mapped_smiles']
        ref_energy = compute_reference_energy(smiles)
        name = name.replace('/', '')  # Remove stereochemistry markers that h5py interprets as path separators
        group = outputfile.create_group(name)
        group.create_dataset('subset', data=[subset], dtype=h5py.string_dtype())
        group.create_dataset('smiles', data=[smiles], dtype=h5py.string_dtype())
        group.create_dataset("atomic_numbers", data=molecules[0].atomic_numbers, dtype=np.int16)
        if max_force is not None:
            force_b3lyp = np.array([vars['DFT TOTAL GRADIENT'] for vars in qcvars_b3lyp])
            force_d3bj = np.array([vars['DISPERSION CORRECTION GRADIENT'] for vars in qcvars_d3bj])
            force = force_b3lyp + force_d3bj
            samples = [i for i in range(len(molecules)) if np.max(np.abs(force[i])) <= max_force]
            molecules = [molecules[i] for i in samples]
            qcvars_b3lyp = [qcvars_b3lyp[i] for i in samples]
            qcvars_d3bj = [qcvars_d3bj[i] for i in samples]
        ds = group.create_dataset('conformations', data=np.array([m.geometry for m in molecules]), dtype=np.float32)
        ds.attrs['units'] = 'bohr'
        ds = group.create_dataset('formation_energy', data=np.array([vars_b3lyp['DFT TOTAL ENERGY']+vars_d3bj['DISPERSION CORRECTION ENERGY']-ref_energy for vars_b3lyp,vars_d3bj in zip(qcvars_b3lyp,qcvars_d3bj)]), dtype=np.float32)
        ds.attrs['units'] = 'hartree'
        for value in config['values']:
            key = value.lower().replace(' ', '_')
            try:
                if key.startswith("dispersion"):
                    qcvars = qcvars_d3bj
                else:
                    qcvars = qcvars_b3lyp
                ds = group.create_dataset(key, data=np.array([vars[value] for vars in qcvars], dtype=np.float32), dtype=np.float32)
                if key in units:
                    ds.attrs['units'] = units[key]
            except KeyError:
                pass
