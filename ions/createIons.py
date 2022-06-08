import numpy as np
import h5py

ions = [('Br', -1), ('Cl', -1), ('F', -1), ('I', -1), ('K', 1), ('Li', 1), ('Na', 1)]

def charge(c):
    if c == 1:
        return '+'
    if c == -1:
        return '-'
    return '%+d' % c

# Create the pairs.

outputfile = h5py.File('ions.hdf5', 'w')
for i in range(len(ions)):
    ion1 = ions[i]
    for j in range(i+1):
        ion2 = ions[j]
        smiles = f'[{ion1[0]}{charge(ion1[1])}:1].[{ion2[0]}{charge(ion2[1])}:2]'
        conformations = []
        for r in range(25, 76):
            conformations.append([[0, 0, 0], [0, 0, 0.01*r]])
        group = outputfile.create_group(smiles)
        group.create_dataset('smiles', data=[smiles], dtype=h5py.string_dtype())
        ds = group.create_dataset('conformations', data=np.array(conformations), dtype=np.float32)
        ds.attrs['units'] = 'nanometers'
