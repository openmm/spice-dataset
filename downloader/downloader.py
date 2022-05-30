from qcportal import FractalClient
from collections import defaultdict
from rdkit import Chem
import numpy as np
import h5py
import yaml

# Reference energies computed with Psi4 1.5 wB97M-D3BJ/def2-TZVPPD.

atom_energy = {('Br', -1): -2574.2451510945853, ('Br', 0): -2574.1167240829964, ('C', -1): -37.91424135791358, ('C', 0): -37.87264507233593,
               ('C', 1): -37.45349214963933, ('Ca', 2): -676.9528465198214, ('Cl', -1): -460.3350243496703, ('Cl', 0): -460.1988762285739,
               ('F', -1): -99.91298732343974, ('F', 0): -99.78611622985483, ('H', 0): -0.4987605100487341, ('I', -1): -297.8813829975981,
               ('I', 0): -297.76228914445625, ('K', 1): -599.8025677513111, ('Li', 1): -7.285254714046546, ('Mg', 2): -199.2688420040449,
               ('N', -1): -54.602291095426494, ('N', 0): -54.62327513368922, ('N', 1): -54.08594142587869, ('Na', 1): -162.11366478783253,
               ('O', -1): -75.17101657391741, ('O', 0): -75.11317840410095, ('O', 1): -74.60241514396725, ('P', 0): -341.3059197024934,
               ('P', 1): -340.9258392474849, ('S', -1): -398.2405387031612, ('S', 0): -398.1599636677874, ('S', 1): -397.7746615977658}

with open('config.yaml') as input:
    config = yaml.safe_load(input.read())
client = FractalClient()
outputfile = h5py.File('SPICE.hdf5', 'w')
for subset in config['subsets']:
    # Download the next subset.

    print('Processing', subset)
    ds = client.get_collection('Dataset', subset)
    all_molecules = ds.get_molecules()
    spec = ds.list_records().iloc[0].to_dict()
    recs = ds.get_records(method=spec['method'], basis=spec['basis'], program=spec['program'], keywords=spec['keywords'])
    recs_by_name = defaultdict(list)
    mols_by_name = defaultdict(list)
    for i in range(len(recs)):
        rec = recs.iloc[i].record
        if rec is not None and rec.status == 'COMPLETE':
            index = recs.index[i]
            name = index[:index.rfind('-')]
            recs_by_name[name].append(rec)
            mols_by_name[name].append(all_molecules.loc[index][0])

    # Add the data to the HDF5 file.

    for name in recs_by_name:
        group_recs = recs_by_name[name]
        molecules = mols_by_name[name]
        qcvars = [r.extras['qcvars'] for r in group_recs]
        smiles = molecules[0].extras['canonical_isomeric_explicit_hydrogen_mapped_smiles']
        rdmol = Chem.MolFromSmiles(smiles, sanitize=False)
        ref_energy = sum(atom_energy[(atom.GetSymbol(), atom.GetFormalCharge())] for atom in rdmol.GetAtoms())
        group = outputfile.create_group(name)
        group.create_dataset('subset', data=[subset], dtype=h5py.string_dtype())
        group.create_dataset('smiles', data=[smiles], dtype=h5py.string_dtype())
        group.create_dataset("atomic_numbers", data=molecules[0].atomic_numbers, dtype=np.int16)
        group.create_dataset('conformations', data=np.array([m.geometry for m in molecules]), dtype=np.float32)
        group.create_dataset('formation_energy', data=np.array([vars['DFT TOTAL ENERGY']-ref_energy for vars in qcvars]), dtype=np.float32)
        for value in config['values']:
            key = value.lower().replace(' ', '_')
            try:
                group.create_dataset(key, data=np.array([vars[value] for vars in qcvars]), dtype=np.float32)
            except KeyError:
                pass
