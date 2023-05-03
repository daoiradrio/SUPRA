import os
import numpy as np
from utils.analyzer import Analyzer
from SUPRAConformer.structure import Structure
from scipy.optimize import linear_sum_assignment



def rmsd_cost(atom1, coord1, atom2, coord2):
    element1 = Structure.get_element(atom1)
    element2 = Structure.get_element(atom2)
    if element1 == element2:
        element_term = 0.0
    else:
        element_term = 100.0
    return (coord1[0] - coord2[0])**2 + (coord1[1] - coord2[1])**2 + (coord1[2] - coord2[2])**2 + element_term


def new_filter_doubles(path: str):
    conformer1 = Structure()
    conformer2 = Structure()
    path = os.path.abspath(path)
    liste = os.listdir(path)
    counter = len(liste)
    for index, file1 in enumerate(liste):
        conformer1.read_xyz(os.path.join(path, file1))
        for file2 in liste[index + 1:]:
            conformer2.read_xyz(os.path.join(path, file2))
            new_coords1, new_coords2 = Analyzer.kabsch(conformer1.coords, conformer2.coords)
            n_atoms = len(conformer1.coords.keys())
            cost = np.zeros((n_atoms, n_atoms))
            for i, (atom1, coords1) in enumerate(zip(conformer1.coords.keys(), new_coords1)):
                for j, (atom2, coords2) in enumerate(zip(conformer2.coords.keys(), new_coords2)):
                    cost[i][j] = rmsd_cost(atom1, coords1, atom2, coords2)
            row, col = linear_sum_assignment(cost)
            final_coords1 = []
            final_coords2 = []
            for i, j in zip(row, col):
                final_coords1.append(new_coords1[i])
                final_coords2.append(new_coords2[j])
            if Analyzer.rmsd(final_coords1, final_coords2) <= 0.1:
                os.remove(os.path.join(path, file1))
                counter -= 1
                break
    print(f"Individual conformers in {path}: {counter}")
    print(f"Total time RMSD: {self.total_time_rmsd}")
    print(f"Total time get_structure: {self.total_time_get_structure}")


a = Analyzer()
a.filter_doubles("/home/baum/SUPRA/tryout/SUPRA_Output/")
#new_filter_doubles("/home/baum/SUPRA/tryout/SUPRA_Output/")
