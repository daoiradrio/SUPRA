from SUPRAConformer.structure import Structure
from utils.helper import get_element, geometric_center
from utils.analyzer import Analyzer
from scipy.spatial.transform import Rotation
import numpy as np



#file1 = "tests/testcases/Alanin.xyz"
#file2 = "tests/testcases/Alanin_different_atom_order.xyz"
file1 = "../kabsch_fail/conformer0.xyz"
file2 = "../kabsch_fail/conformer1.xyz"

struc1 = Structure(file1)
struc2 = Structure(file2)

analyzer = Analyzer()
print(analyzer._rmsd_tight(struc1.coords, struc1.bond_partners, struc2.coords, struc2.bond_partners))

new_coords1, new_coords2 = analyzer.SCHK(struc1.coords, struc2.coords)
print(analyzer._calc_rmsd(new_coords1, new_coords2))

"""
benzene = Structure("tests/testcases/Benzol.xyz")
angle = 180.0
geomcenter = geometric_center(benzene.coords)
v1 = benzene.coords["C0"] - geomcenter
v1 = v1 / np.linalg.norm(v1)
v2 = benzene.coords["C1"] - geomcenter
v2 = v2 / np.linalg.norm(v2)
norm_vec = np.cross(v1, v2)
norm_vec = norm_vec / np.linalg.norm(norm_vec)
norm_vec = geomcenter + norm_vec
R = Rotation.from_rotvec(np.deg2rad(angle) * norm_vec)
for atom, coord in benzene.coords.items():
    #print()
    #print(coord)
    coord = R.apply(coord)
    x, y, z = coord
    #print(coord)
    print(f"{get_element(atom)}\t{x:18.15f}\t{y:18.15f}\t{z:18.15f}")
"""
