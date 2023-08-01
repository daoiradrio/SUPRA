import numpy as np

from utils.analyzer import Analyzer
from SUPRAConformer.structure import Structure



analyzer = Analyzer()
mol1 = Structure("tests/testcases/Alanin.xyz")
mol2 = Structure("tests/testcases/Alanin_different_atom_order.xyz")

coords1 = np.array(list(mol1.coords.values()))
coords2 = np.array(list(mol2.coords.values()))

new_coords1, new_coords2 = analyzer.SCHK(mol1.coords, mol2.coords)

print(analyzer._calc_rmsd(new_coords1, new_coords2))