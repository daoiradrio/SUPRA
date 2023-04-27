import os
import numpy as np
from helper import inertia_tensor
from SUPRAConformer.structure import Structure



path = "/home/baum/SUPRA/inputfiles/"
xyz_file = "Alanin.xyz"
xyz_path = os.path.join(path, xyz_file)
mol = Structure()
mol.read_xyz(xyz_path)
I = inertia_tensor(mol.coords)
eig_vals, eig_vecs = np.linalg.eig(I)
e1, e2, e3 = eig_vals
tol = min(e1, e2, e3)
delta = min(abs(e1-e2), abs(e1-e3), abs(e2-e3))
print(eig_vals)
print(delta)
print(tol)
if delta <= tol:
    print("Symmetrical")