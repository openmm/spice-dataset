import openmm.app as app
from rdkit import Chem
from concurrent.futures import ThreadPoolExecutor

from rdkit import RDLogger
RDLogger.DisableLog('rdApp.*')

elements = set(['B', 'Br', 'C', 'Cl', 'F', 'H', 'I', 'N', 'O', 'P', 'S', 'Si'])

def processLigand(smiles, resid):
    try:
        mol = Chem.MolFromSmiles(smiles)
        mol = Chem.AddHs(mol)
        from rdkit.Chem.MolStandardize import rdMolStandardize
        enumerator = rdMolStandardize.TautomerEnumerator()
        enumerator.SetMaxTautomers(2)
        results = []
        for variant in enumerator.Enumerate(mol):
            results.append(Chem.MolToSmiles(variant))
        return resid, results
    except:
        return None

futures = []
with ThreadPoolExecutor() as executor:
    for line in open('Components-smiles-oe.smi'):
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
        if mol.GetNumAtoms() < 3 or mol.GetNumAtoms() > 100:
            continue
        futures.append(executor.submit(processLigand, smiles, resid))

allSmiles = set()
with open('ligandExpoVariants.txt', 'w') as output:
    for future in futures:
        resid, variants = future.result()
        for i, smiles in enumerate(variants):
            if smiles not in allSmiles:
                print(f'{resid} variant {i+1}\t{smiles}', file=output)
                allSmiles.add(smiles)