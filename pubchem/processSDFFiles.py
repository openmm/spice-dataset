from rdkit import Chem
import os

elements = set(['B', 'Br', 'C', 'Ca', 'Cl', 'F', 'H', 'I', 'K', 'Li', 'Mg', 'N', 'Na', "O", 'P', 'S', 'Si'])

for filename in os.listdir('sources'):
    if filename.endswith('.sdf'):
        print('Processing', filename)
        inputFile = os.path.join('sources', filename)
        outputFile = inputFile[:-4]+'.txt'
        mols = Chem.SDMolSupplier(inputFile)
        total = 0
        included = 0
        with open(outputFile, 'w') as output:
            for mol in mols:
                total += 1
                if mol is None or mol.GetNumAtoms() < 3:
                    continue
                if any(a.GetSymbol() not in elements for a in mol.GetAtoms()):
                    continue
                try:
                    mol2 = Chem.AddHs(mol)
                except:
                    continue
                if mol2.GetNumAtoms() > 50:
                    continue
                if any(a.GetNumRadicalElectrons() != 0 or a.GetIsotope() != 0 for a in mol2.GetAtoms()):
                    continue
                smiles = Chem.MolToSmiles(mol)
                if '.' in smiles:
                    continue
                molid = mol.GetProp('PUBCHEM_SUBSTANCE_ID')
                output.write(f'{molid}\t{smiles}\n')
                included += 1
        print('  Included', included, 'of', total)
