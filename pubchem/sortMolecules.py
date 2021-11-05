from rdkit import Chem, DataStructs
from rdkit.Chem import AllChem
import numpy as np
import os

class Molecule(object):
    def __init__(self, sid, smiles):
        self.sid = sid
        self.smiles = smiles
        mol = Chem.MolFromSmiles(smiles)
        if any(a.GetSymbol() in ('B', 'Si') for a in mol.GetAtoms()):
            raise ValueError('Skipping molecules with B or Si')
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

    count = 0
    while len(fp) > 0:
        i = np.argmin(maxSimilarity)
        if count%1000 == 0:
            print(count, len(fp), maxSimilarity[i])
        m = mols[i]
        del mols[i]
        del fp[i]
        maxSimilarity = np.delete(maxSimilarity, i)
        output.write(f'{m.sid}\t{m.smiles}')
        similarity = DataStructs.BulkTanimotoSimilarity(m.fingerprint, fp)
        maxSimilarity = np.maximum(similarity, maxSimilarity)
        if count%1000 == 0:
            # Remove any molecule that is too similar to one we have already added.  Including them
            # wouldn't significantly expand our coverage of chemical space.

            ones = (maxSimilarity >= 0.9).nonzero()[0]
            for j in reversed(ones):
                del mols[j]
                del fp[j]
            maxSimilarity = np.delete(maxSimilarity, ones)
        count += 1