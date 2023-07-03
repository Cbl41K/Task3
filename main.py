import pandas as pd
from rdkit import Chem
from rdkit.Chem import Descriptors

df = pd.read_csv('SolubilityDataset.csv')

data = []


for _, row in df.iterrows():
    smiles = row['SMILES']
    mol = Chem.MolFromSmiles(smiles)

    if mol is not None:
        solubility = row['Solubility']
        mol_log_p = Descriptors.MolLogP(mol)
        mol_mr = Descriptors.MolMR(mol)
        heavy_atom_count = Descriptors.HeavyAtomCount(mol)
        num_rotatable_bonds = Descriptors.NumRotatableBonds(mol)
        num_valence_electrons = Descriptors.NumValenceElectrons(mol)
        num_aromatic_rings = Descriptors.NumAromaticRings(mol)
        tpsa = Descriptors.TPSA(mol)
        labute_asa = Descriptors.LabuteASA(mol)
        bertz_ct = Descriptors.BertzCT(mol)

        data.append([smiles, solubility, mol_log_p, mol_mr, heavy_atom_count, num_rotatable_bonds,
                     num_valence_electrons, num_aromatic_rings, tpsa, labute_asa, bertz_ct])

new_df = pd.DataFrame(data, columns=['SMILES', 'Solubility', 'MolLogP', 'MolMR', 'HeavyAtomCount',
                                     'NumRotatableBonds', 'NumValenceElectrons', 'NumAromaticRings',
                                     'TPSA', 'LabuteASA', 'BertzCT'])


new_df.to_csv('DataFrame_new.csv', index=False)