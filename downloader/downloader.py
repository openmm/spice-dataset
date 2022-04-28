from qcportal import FractalClient
from collections import defaultdict
import numpy as np
import h5py
import yaml

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
        if rec.status == 'COMPLETE':
            index = recs.index[i]
            name = index[:index.rfind('-')]
            recs_by_name[name].append(rec)
            mols_by_name[name].append(all_molecules.loc[index][0])

    # Add the data to the HDF5 file.

    for name in recs_by_name:
        group_recs = recs_by_name[name]
        molecules = mols_by_name[name]
        qcvars = [r.extras['qcvars'] for r in group_recs]
        group = outputfile.create_group(name)
        group.create_dataset('subset', data=[subset], dtype=h5py.string_dtype())
        group.create_dataset('smiles', data=[molecules[0].extras['canonical_isomeric_explicit_hydrogen_mapped_smiles']], dtype=h5py.string_dtype())
        group.create_dataset("atomic_numbers", data=molecules[0].atomic_numbers, dtype=np.int16)
        conformations = group.create_dataset('conformations', data=np.array([m.geometry for m in molecules]), dtype=np.float32)
        for value in config['values']:
            key = value.lower().replace(' ', '_')
            try:
                group.create_dataset(key, data=np.array([vars[value] for vars in qcvars]), dtype=np.float32)
            except KeyError:
                pass
