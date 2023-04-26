import os
import numpy as np
from scipy.optimize import linear_sum_assignment
from SUPRAConformer.structure import Structure
from utils.analyzer import Analyzer



def kabsch(coords1, coords2):
    # KABSCH ALOGIRTHMUS
    # Minimiere RMSD zwischen zwei Koordinatensätzen
    
    # Falls nötig Datentypen anpassen
    if type(coords1) == dict:
        coords1 = list(coords1.values())
        coords1 = np.array(coords1)
    if type(coords2) is dict:
        coords2 = list(coords2.values())
        coords2 = np.array(coords2)

    # Zentriere beide Sets von Koordinaten über jeweilige geometrische Zentren
    center1 = np.mean(coords1, axis=0)
    center2 = np.mean(coords2, axis=0)
    coords1 -= center1
    coords2 -= center2

    # Kovarianzmatrix berechnen
    H = np.matmul(coords1.T, coords2)

    # Singulärwertzerlegung der Kovarianzmatrix berechnen
    U, S, Vt = np.linalg.svd(H)

    # Matrix zur Berechnung der Rotationsmatrix in Abhängigkeit der Determinante bestimmen
    det = np.linalg.det(np.matmul(Vt.T, U.T))
    if det >= 0:
        det = 1.0
    else:
        det = -1.0
    matrix = np.array([[1, 0, 0], [0, 1, 0], [0, 0, det]])

    # Rotationsmatrix berechnen
    R = np.matmul(np.matmul(Vt.T, matrix), U.T)

    # anwenden der Rotationsmatrix auf Koordinatenset 2 um beide Sets möglichst zur Deckung zu bringen
    for index in range(len(coords2)):
        coords2[index] = np.matmul(coords2[index], R)
    
    return (coords1, coords2)



def write_xyz(coords: dict):
    with open("out.xyz", "w") as xyz_file:
        n_atoms = len(coords.keys())
        print(n_atoms, file=xyz_file, end="\n\n")
        for atom, coord in coords.items():
            element = Structure.get_element(atom)
            x = coord[0]
            y = coord[1]
            z = coord[2]
            print(f"{element}\t{x}\t{y}\t{z}", file=xyz_file)



def rmsd_cost(atom1, coord1, atom2, coord2):
    element1 = Structure.get_element(atom1)
    element2 = Structure.get_element(atom2)
    if element1 == element2:
        element_term = 0.0
    else:
        element_term = 100.0
    return (coord1[0] - coord2[0])**2 + (coord1[1] - coord2[1])**2 + (coord1[2] - coord2[2])**2 + element_term



folder_testcases = "/home/dario/SUPRA/tests/testcases/"
file1 = "Tyrosin.xyz"
file2 = "Tyrosin_phenyl_rotated_180_deg.xyz"
path1 = os.path.join(folder_testcases, file1)
path2 = os.path.join(folder_testcases, file2)
mol1 = Structure(path1)
mol2 = Structure(path2)
a = Analyzer()

new_coords1, new_coords2 = Analyzer.kabsch(mol1.coords, mol2.coords)
n_atoms = len(mol1.coords.keys())
cost = np.zeros((n_atoms, n_atoms))
for i, (atom1, coords1) in enumerate(zip(mol1.coords.keys(), new_coords1)):
    for j, (atom2, coords2) in enumerate(zip(mol2.coords.keys(), new_coords2)):
        cost[i][j] = rmsd_cost(atom1, coords1, atom2, coords2)
row, col = linear_sum_assignment(cost)
final_coords1 = []
final_coords2 = []
for i, j in zip(row, col):
    final_coords1.append(new_coords1[i])
    final_coords2.append(new_coords2[j])
print(Analyzer.rmsd(final_coords1, final_coords2))

new_coords1, new_coords2 = a.match(
        mol1.coords, mol1.bond_partners,
        mol2.coords, mol2.bond_partners
)
print(Analyzer.kabsch_and_rmsd(new_coords1, new_coords2))
