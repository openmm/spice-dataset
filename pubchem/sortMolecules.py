from rdkit import Chem, DataStructs
from rdkit.Chem import AllChem
import numpy as np
import os

class Molecule(object):
    def __init__(self, sid, smiles):
        self.sid = sid
        self.smiles = smiles
        mol = Chem.MolFromSmiles(smiles)
        mol = Chem.AddHs(mol)
        self.fingerprint = AllChem.GetMorganFingerprintAsBitVect(mol, 2, 1024)

mols = []

for filename in os.listdir('sources'):
    if filename.endswith('.txt'):
        print('Loading', filename)
        inputFile = os.path.join('sources', filename)
        for line in open(inputFile):
            sid, smiles = line.split('\t')
            try:
                mols.append(Molecule(sid, smiles))
            except:
                pass

print('Sorting')
with open('sorted.txt', 'w') as output:
    # Start with the first molecule.

    output.write(f'{mols[0].sid}\t{mols[0].smiles}')
    fp = [m.fingerprint for m in mols[1:]]
    maxSimilarity = DataStructs.BulkTanimotoSimilarity(mols[0].fingerprint, fp)
    del mols[0]

    # Repeatedly add the molecule least similar to anything so far.

    while len(fp) > 0:
        i = np.argmin(maxSimilarity)
        m = mols[i]
        del mols[i]
        del fp[i]
        maxSimilarity = np.delete(maxSimilarity, i)
        output.write(f'{m.sid}\t{m.smiles}')
        similarity = DataStructs.BulkTanimotoSimilarity(m.fingerprint, fp)
        maxSimilarity = np.maximum(similarity, maxSimilarity)
