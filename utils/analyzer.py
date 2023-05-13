import os

import time

import numpy as np
import multiprocessing as mp

from utils.helper import get_element, valences
from SUPRAConformer.structure import Structure
from typing import Union
from queue import Queue
from scipy.optimize import linear_sum_assignment



class Analyzer:

    def __init__(self):
        self.connectivity = {}
        self.ignored_terminal_group_atoms = []
        self.ignored_methyl_group_atoms = []


    def compare_structure_sets(self, path1: str, path2: str):
        Molecule1 = Structure()
        Molecule2 = Structure()
        list1 = os.listdir(path1)
        list2 = os.listdir(path2)
        doubles = 0
        for i, conformer1 in enumerate(list1):
            print("Progress {:.1f}%".format(i/len(list1) * 100))
            Molecule1.get_structure(path1 + conformer1)
            for conformer2 in list2:
                Molecule2.get_structure(path2 + conformer2)
                if self.doubles(Molecule1, Molecule2):
                    doubles += 1
                    break
        print("Progress 100.0%")
        print()
        print("Number of structures in " + path1 + " (path 1): " + str(len(list1)))
        print("Number of structures in " + path2 + " (path 2): " + str(len(list2)))
        print("Number of structures of path 1 in path 2: " + str(doubles))

    
    # TODO: ÜBERARBEITEN, AUF NEUESTEN STAND BRINGEN
    #def filter_doubles_parallel(self, path: str, cpu: int):
    #    if cpu > mp.cpu_count():
    #        cpu = mp.cpu_count()
    #        print("Exceeds maxmimum number of CPU, "
    #              "calculation will be running with maximum number of " + str(cpu) + " CPU.")
    #    amount = len(os.listdir(path))
    #    processes = mp.Pool(cpu)
    #    results = processes.starmap(self.filter,
    #                                [(path, int(((n-1)/4)*amount), int((n/4)*amount)) for n in range(1,cpu+1)])
    #    processes.close()
    #    individuals = 0
    #    with open(path[:-1]+".xyz", "w") as outfile:
    #        for result in results:
    #            for conformer in result:
    #                individuals += 1
    #                print(conformer, file=outfile)
    #    print("Individual conformers in " + path + ": " + str(individuals))
    #
    #
    #def filter(self, path: str, min_index: int, max_index: int) -> list:
    #    conformer1 = Structure()
    #    conformer2 = Structure()
    #    conformer_list = os.listdir(path)
    #    conformers = [conformer for conformer in conformer_list[min_index : max_index]]
    #    for iter, file1 in enumerate(conformer_list[min_index : max_index]):
    #        iter = min_index + iter + 1
    #        conformer1.get_structure(path + file1)
    #        for file2 in conformer_list[iter : ]:
    #            conformer2.get_structure(path + file2)
    #            if self.doubles(conformer1, conformer2):
    #                conformers.remove(file1)
    #                os.remove(path + file1)
    #                break
    #    return conformers

    
    def remove_doubles(self, path: str, rmsd_threshold: float=None, ignore: str=None) -> None:
        print("Performing removal of duplicate structures...")

        conformer1 = Structure()
        conformer2 = Structure()
        path = os.path.abspath(path)
        conformers = os.listdir(path)
        counter = len(conformers)

        if not rmsd_threshold:
            rmsd_threshold = 0.1

        if ignore=="methyl":
            conformer1.get_structure(os.path.join(path, conformers[0]))
            self.connectivity = conformer1.bond_partners
            self.set_methyl_group_atoms()
            atoms = [atom for atom in conformer1.coords.keys() if not atom in self.ignored_methyl_group_atoms]
        elif ignore=="all":
            conformer1.get_structure(os.path.join(path, conformers[0]))
            self.connectivity = conformer1.bond_partners
            self.set_terminal_group_atoms()
            atoms = [atom for atom in conformer1.coords.keys() if not atom in self.ignored_terminal_group_atoms]
        else:
            conformer1.read_xyz(os.path.join(path, conformers[0]))
            atoms = [atom for atom in conformer1.coords.keys()]

        n_atoms = len(atoms)
        cost = np.zeros((n_atoms, n_atoms))

        for index, file1 in enumerate(conformers):
            conformer1.read_xyz(os.path.join(path, file1))
            if ignore == "methyl":
                for atom in self.ignored_methyl_group_atoms:
                    del conformer1.coords[atom]
            elif ignore == "all": 
                for atom in self.ignored_methyl_group_atoms:
                    del conformer1.coords[atom]
            for file2 in conformers[index + 1:]:
                conformer2.read_xyz(os.path.join(path, file2))
                if ignore == "methyl": 
                    for atom in self.ignored_methyl_group_atoms:
                        del conformer2.coords[atom]
                elif ignore == "all": 
                    for atom in self.ignored_methyl_group_atoms:
                        del conformer2.coords[atom]
                kabsch_coords1, kabsch_coords2 = self.kabsch(conformer1.coords, conformer2.coords)
                for i in range(n_atoms):
                    for j in range(i+1):
                        element_i = get_element(atoms[i])
                        element_j = get_element(atoms[j])
                        if element_i == element_j:
                            element_term = 0.0 
                        else:
                            element_term = 100.0
                        diff_vec = kabsch_coords1[i] - kabsch_coords2[j]
                        cost_value = np.dot(diff_vec, diff_vec) + element_term
                        cost[i][j] = cost_value
                        cost[j][i] = cost_value
                row, col = linear_sum_assignment(cost)
                print("Anfang Hungarian")
                print()
                print(kabsch_coords1)
                print()
                print(row)
                print()
                print(col)
                print()
                print("Ende Hungarian")
                print()
                if self.rmsd(kabsch_coords1[row], kabsch_coords2[col]) <= rmsd_threshold:
                    #os.remove(os.path.join(path, file1))
                    counter -= 1
                    break
                    
        print(f"Individual conformers in {path}: {counter}")
    

    def set_methyl_group_atoms(self) -> None:
        for atom, bond_partners in self.connectivity.items():
            if get_element(atom) == "C":
                terminal_count = 0
                bond_partners = [get_element(bond_partner) for bond_partner in bond_partners]
                terminal_count += bond_partners.count("H")
                terminal_count += bond_partners.count("F")
                terminal_count += bond_partners.count("Cl")
                terminal_count += bond_partners.count("Br")
                terminal_count += bond_partners.count("I")
                if terminal_count == valences[get_element(atom)]-1:
                    sorted_bond_partners = sorted(self.connectivity[atom], key=lambda x: valences[get_element(x)], reverse=True)
                    self.ignored_methyl_group_atoms.append(atom)
                    self.ignored_methyl_group_atoms.append(sorted_bond_partners[1])
                    self.ignored_methyl_group_atoms.append(sorted_bond_partners[2])
                    self.ignored_methyl_group_atoms.append(sorted_bond_partners[3])


    def set_terminal_group_atoms(self) -> None:
        for atom, bond_partners in self.connectivity.items():
            if not get_element(atom) in ["H", "F", "Cl", "Br", "I"]:
                terminal_count = 0
                bond_partners = [get_element(bond_partner) for bond_partner in bond_partners]
                terminal_count += bond_partners.count("H")
                terminal_count += bond_partners.count("F")
                terminal_count += bond_partners.count("Cl")
                terminal_count += bond_partners.count("Br")
                terminal_count += bond_partners.count("I")
                if terminal_count == valences[get_element(atom)]-1:
                    sorted_bond_partners = sorted(self.connectivity[atom], key=lambda x: valences[get_element(x)], reverse=True)
                    self.ignored_terminal_group_atoms.append(atom)
                    for terminal_atom in sorted_bond_partners[1:]:
                        self.ignored_terminal_group_atoms.append(terminal_atom)


    def kabsch(self, coords1: Union[dict, list, np.array], coords2: Union[dict, list, np.array]) -> tuple:
        if type(coords1) == dict:
            coords1 = list(coords1.values())
            coords1 = np.array(coords1)
        if type(coords2) is dict:
            coords2 = list(coords2.values())
            coords2 = np.array(coords2)
        center1 = np.mean(coords1, axis=0)
        center2 = np.mean(coords2, axis=0)
        coords1 -= center1
        coords2 -= center2

        # Kovarianzmatrix berechnen
        H = np.matmul(coords1.T, coords2)
        #print("H:")
        #print(H)
        #print()

        # Singulärwertzerlegung der Kovarianzmatrix berechnen
        U, S, Vt = np.linalg.svd(H)
        #print("U:")
        #print(U)
        #print()
        #print("S:")
        #print(S)
        #print()
        #print("Vt:")
        #print(Vt)
        #print()

        # Matrix zur Berechnung der Rotationsmatrix in Abhängigkeit der Determinante bestimmen
        det = np.linalg.det(np.matmul(Vt.T, U.T))
        #print("det:")
        #print(det)
        #print()
        if det >= 0:
            det = 1.0
        else:
            det = -1.0
        matrix = np.array([[1, 0, 0], [0, 1, 0], [0, 0, det]])

        # Rotationsmatrix berechnen
        R = np.matmul(np.matmul(Vt.T, matrix), U.T)

        #print("R:")
        #print(R)
        #print()

        #print("coords2 vorher:")
        #print(coords2)
        #print()

        # anwenden der Rotationsmatrix auf Koordinatenset 2 um beide Sets möglichst zur Deckung zu bringen
        for i, _ in enumerate(coords2):
            coords2[i] = np.matmul(coords2[i], R)

        #print("coords2 nachher:")
        #print(coords2)
        #print()

        return (coords1, coords2)
    

    def rmsd(self, coords1: Union[dict, list, np.array], coords2: Union[dict, list, np.array]) -> float:
        if type(coords1) == dict:
            coords1 = list(coords1.values())
        if type(coords2) is dict:
            coords2 = list(coords2.values())
        n = len(coords1)
        delta_sum = 0.0
        for i in range(n):
            delta_sum += (coords1[i][0] - coords2[i][0])**2 + \
                         (coords1[i][1] - coords2[i][1])**2 + \
                         (coords1[i][2] - coords2[i][2])**2
        return np.sqrt(1.0/float(n) * delta_sum)


    # überprüfen, ob zwei Strukturen identisch (Doubles) oder verschieden sind
    #def rmsd(self, coords1: Union[Structure, dict, list], coords2: Union[Structure, dict, list]) -> float:
    @staticmethod
    def kabsch_and_rmsd(coords1: Union[Structure, dict, list], coords2: Union[Structure, dict, list]) -> float:
        # ggf. umwandeln des Koordinatendateityps von Dictionary zu Liste für beide Koordinatensets
        if type(coords1) == Structure:
            coords1 = list(coords1.coords.values())
        elif type(coords1) == dict:
            coords1 = list(coords1.values())

        if type(coords2) == Structure:
            coords2 = list(coords2.coords.values())
        elif type(coords2) is dict:
            coords2 = list(coords2.values())

        # Überprüfen ob beide Koordinatensets dieselbe Anzahl an Punkten enthalten
        if len(coords1) != len(coords2):
            print("Number of atoms are not the same in the given structures.")
            raise ValueError
        else:
            coords1 = np.array(coords1)
            coords2 = np.array(coords2)

        # 1. Minimiere RMSD zwischen self.coords und coords mithilfe von Kabsch-Algorithmus

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

        # 2. Berechne RMSD

        # Anzahl der Punkte N, es gilt len(coords1)=len(coords2) (s.o.)
        N = len(coords1)
        sum = 0.0
        for i in range(N):
            sum += (coords1[i][0] - coords2[i][0])**2 + \
                   (coords1[i][1] - coords2[i][1])**2 + \
                   (coords1[i][2] - coords2[i][2])**2
        rmsd = np.sqrt(1.0/float(N) * sum)

        return rmsd
