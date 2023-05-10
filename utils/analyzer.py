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


    def filter_doubles_parallel(self, path: str, cpu: int):
        if cpu > mp.cpu_count():
            cpu = mp.cpu_count()
            print("Exceeds maxmimum number of CPU, "
                  "calculation will be running with maximum number of " + str(cpu) + " CPU.")
        amount = len(os.listdir(path))
        processes = mp.Pool(cpu)
        results = processes.starmap(self.filter,
                                    [(path, int(((n-1)/4)*amount), int((n/4)*amount)) for n in range(1,cpu+1)])
        processes.close()
        individuals = 0
        with open(path[:-1]+".xyz", "w") as outfile:
            for result in results:
                for conformer in result:
                    individuals += 1
                    print(conformer, file=outfile)
        print("Individual conformers in " + path + ": " + str(individuals))


    def filter(self, path: str, min_index: int, max_index: int) -> list:
        conformer1 = Structure()
        conformer2 = Structure()
        conformer_list = os.listdir(path)
        conformers = [conformer for conformer in conformer_list[min_index : max_index]]
        for iter, file1 in enumerate(conformer_list[min_index : max_index]):
            iter = min_index + iter + 1
            conformer1.get_structure(path + file1)
            for file2 in conformer_list[iter : ]:
                conformer2.get_structure(path + file2)
                if self.doubles(conformer1, conformer2):
                    conformers.remove(file1)
                    os.remove(path + file1)
                    break
        return conformers


    def filter_doubles_multiple(self, path: str, molecule: str):
        conformer1 = Structure()
        conformer2 = Structure()
        for folder in os.listdir(path):
            if not molecule in folder:
                continue
            extension = path + folder + "/"
            liste = os.listdir(extension)
            counter = len(liste)
            conformers = []
            for index, file1 in enumerate(liste):
                print(folder + ": " + "{:.2f}".format((index+1)*(100/len(liste))) + "%")
                #conformer1.get_structure(extension + file1)
                conformer1.get_structure(extension + file1)
                flag = False
                for file2 in liste[index + 1:]:
                    conformer2.get_structure(extension + file2)
                    if self.doubles(conformer1, conformer2):
                        counter -= 1
                        #os.remove(extension + file1)
                        flag = True
                        break
                if not flag:
                    conformers.append(file1)
            with open(path+folder+".xyz", "w") as outfile:
                for conformer in conformers:
                    print(conformer, file=outfile)
            #print(folder + ": " + str(counter))
    

    # No. 1 are reference coordinates and structure
    # No. 2 coordinates will be reordered to match No. 1
    def match(self, coords1: dict, structure1: dict, coords2: dict, structure2: dict) -> tuple:
        # initialization
        pairs = {}
        matched_coords1 = []
        matched_coords2 = []
        list1 = list(coords1.keys())
        list2 = list(coords2.keys())
        spheres1 = {atom: self.spheres(structure1, atom) for atom in list1}
        spheres2 = {atom: self.spheres(structure2, atom) for atom in list2}

        # find pairs of atoms from No. 1 and 2 that are equivalent in terms of element and coordination spheres
        # for atoms where this is already unique store the respective coordinates
        for atom1 in list1:
            eq = []
            for atom2 in list2:
                # compare element symbol
                if get_element(atom1) == get_element(atom2):
                    # compare coordination spheres (which elements are in which coordination sphere)
                    if spheres1[atom1] == spheres2[atom2]:
                        #pairs[atom1].append(atom2)
                        eq.append(atom2)
            # if several atoms from No. 2 are equivalent to a respective atom of No. 1 in terms of element symbol and
            # coordination spheres store these candidates
            if len(eq) > 1:
                pairs[atom1] = eq
            # if an atom pair is already unique store the respective coordinates
            else:
                matched_coords1.append(coords1[atom1])
                matched_coords2.append(coords2[eq[0]])

        # find correct atom pairs from No. 1 and No. 2 via matching of distances to other atoms
        for atom, eq_atoms in pairs.items():
            d1 = np.zeros(shape=(len(matched_coords1), len(eq_atoms)))
            d2 = np.zeros(shape=(len(matched_coords2), len(eq_atoms)))
            # store distances Pos.(atom No. 1)-Pos.(matched atom i from No. 1) in matrix d1 which is
            # of shape (number of already matched coordinates)x(number of candidate atoms) which means
            # every row stands for a reference atom from No. 1, every column contains the distance of the
            # yet unmatched atom to the respective reference atom
            for i, comp_atom in enumerate(matched_coords1):
                v11 = np.array(coords1[atom])
                v12 = np.array(comp_atom)
                d1[i] = np.linalg.norm(v11 - v12)
                # store reference distances Pos.(atom No. 2)-Pos.(matched atom i from No. 2) in matrix d2 which is
                # of shape (number of already matched coordinates)x(number of candidate atoms) which means
                # every row stands for a reference atom from No. 2, every column contains the distance of the
                # candidate atom to the respective reference atom
                for j, eq_atom in enumerate(eq_atoms):
                    v21 = np.array(coords2[eq_atom])
                    v22 = np.array(matched_coords2[i])
                    d2[i][j] = np.linalg.norm(v21 - v22)
            # evaluate deviations of distances, after transposing every row stands for a candidate atom and
            # contains columnwise the respective deviations
            d2 = np.transpose(d2 - d1)
            min_delta = 1000000.0
            eq = ""
            # evaluate candidate atom with smallest deviation (compute sum over columns)
            for i, d in enumerate(d2):
                delta = np.linalg.norm(d)
                if delta < min_delta:
                    min_delta = delta
                    eq = pairs[atom][i]
            # store coordinates of atom from No. 1 and coordinates of atom from No. 2 with smallest deviation
            matched_coords1.append(coords1[atom])
            matched_coords2.append(coords2[eq])
            # remove candidate from candidate lists
            for candidates in pairs.values():
                if eq in candidates:
                    candidates.remove(eq)

        return matched_coords1, matched_coords2


    # TODO: this may be more efficient than current version? test on large systems!
    def match_old(self, molecule1: Structure, molecule2: Structure) -> tuple:
        list1 = list(molecule1.coords.keys())
        #list1 = sorted(list1, key=Helper.get_element)
        list2 = list(molecule2.coords.keys())
        #list2 = sorted(list2, key=Helper.get_element)
        pairs = {atom: [] for atom in list1}
        for atom1 in list1:
            element1 = get_element(atom1)
            for atom2 in list2:
                element2 = get_element(atom2)
                if element1 == element2:
                    #ANPASSEN, SODASS DAS NICHT FÜR JEDES ATOM MEHRMALS GEMACHT WERDEN MUSS
                    if self.spheres(molecule1.bond_partners, atom1) == self.spheres(molecule2.bond_partners, atom2):
                        pairs[atom1].append(atom2)
                #WENN LISTEN VORHER SORTIERT WURDEN (S.O.)
                #else:
                #   break
        coords1 = []
        coords2 = []
        for atom, eq_atoms in pairs.items():
            if len(eq_atoms) > 1:
                vec11 = np.array(molecule1.coords[atom])
                candidates = [[candidate, np.array(molecule2.coords[candidate]), 0] for candidate in eq_atoms]
                for comp_atom, eq_comp_atoms in pairs.items():
                    if len(eq_comp_atoms) == 1:
                        vec12 = np.array(molecule1.coords[comp_atom])
                        distance1 = np.linalg.norm(vec11 - vec12)
                        vec22 = np.array(molecule2.coords[eq_comp_atoms[0]])
                        for candidate in candidates:
                            distance2 = np.linalg.norm(candidate[1] - vec22)
                            candidate[2] += np.abs(distance1 - distance2)
                candidates = sorted(candidates, key=lambda candidate: candidate[2])
                eq_atom = candidates[0][0]
                pairs[atom] = [eq_atom]
                for uneq_atom in pairs:
                    if not uneq_atom == atom and eq_atom in pairs[uneq_atom]:
                        pairs[uneq_atom].remove(eq_atom)
            coords1.append(molecule1.coords[atom])
            coords2.append(molecule2.coords[pairs[atom][0]])
        return (coords1, coords2)


    def spheres(self, structure: dict, start: str) -> list:
        atoms = Queue()
        for atom in structure[start]:
            atoms.put(atom)

        memory = [start]
        spheres = []

        sphere_counter = len(structure[start])
        next_sphere_counter = 0
        iter = 0

        while not atoms.empty():
            spheres.append([])
            while sphere_counter:
                atom = atoms.get()
                spheres[iter].append(get_element(atom))
                sphere_counter -= 1
                for neighbor in structure[atom]:
                    if not neighbor in memory:
                        memory.append(neighbor)
                        atoms.put(neighbor)
                        next_sphere_counter += 1
            spheres[iter] = sorted(spheres[iter])
            iter += 1
            sphere_counter = next_sphere_counter
            next_sphere_counter = 0

        return spheres


    def spheres_old(self, structure: dict, start: str) -> list:
        atoms = Queue()
        atoms.put(start)
        memory = []
        memory.append(start)
        sphere = []
        spheres = []
        sphere_counter = len(structure[start])
        next_sphere_counter = 0
        while not atoms.empty():# and not sphere_counter == 0: <-- WOZU DAS HIER IN ALTER VERSION??:
            atom = atoms.get()
            for neighbor in structure[atom]:
                if neighbor not in memory:
                    memory.append(neighbor)
                    atoms.put(neighbor)
                    sphere.append(get_element(neighbor))
                    sphere_counter -= 1
                    for atom in structure[neighbor]:
                        if not atom in memory:
                            next_sphere_counter += 1
                #DAS HIER SOLLTE NICHT INNERHALB DES OBIGEN LOOPS PASSIEREN, DA ALLE NEIGHBORS ZWANGSLÄUFIG ZU EINER
                #SPHERE GEHÖREN
                if sphere_counter == 0:
                    spheres.append(sorted(sphere))
                    sphere = []
                    sphere_counter = next_sphere_counter
                    next_sphere_counter = 0
        return spheres

    
    def remove_doubles(self, path: str, rmsd_threshold: float=0.1, loose: bool=False) -> None:
        print("Performing removal of duplicate structures...")

        conformer1 = Structure()
        conformer2 = Structure()
        path = os.path.abspath(path)
        conformers = os.listdir(path)
        counter = len(conformers)
        if loose:
            conformer1.get_structure(os.path.join(path, conformers[0]))
            self.connectivity = conformer1.bond_partners
            self.set_terminal_group_atoms()
        else:
            conformer1.read_xyz(os.path.join(path, conformers[0]))
        atoms = [atom for atom in conformer1.coords.keys()]
        n_atoms = len(atoms)
        cost = np.zeros((n_atoms, n_atoms))

        for index, file1 in enumerate(conformers):
            conformer1.read_xyz(os.path.join(path, file1))
            for file2 in conformers[index + 1:]:
                conformer2.read_xyz(os.path.join(path, file2))
                # VERBESSERUNGSPOTENTIAL BZGL. IMPLEMENTIERUNG AB HIER??
                kabsch_coords1, kabsch_coords2 = self.kabsch(conformer1.coords, conformer2.coords)
                if loose:
                    for atom in self.ignored_terminal_group_atoms:
                        del kabsch_coords1[atom]
                        del kabsch_coords2[atom]
                for i in range(n_atoms-1):
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
                if self.rmsd(kabsch_coords1[row], kabsch_coords2[col]) <= rmsd_threshold:
                    os.remove(os.path.join(path, file1))
                    counter -= 1
                    break
                    
        print(f"Individual conformers in {path}: {counter}")


    def set_terminal_group_atoms(self) -> None:
        for atom, bond_partners in self.connectivity.items():
            terminal_count = [get_element(bond_partner) for bond_partner in bond_partners].count(["H", "F", "Cl", "Br", "I"])
            if terminal_count == 3:
                sorted_bond_partners = sorted(self.connectivity, key=lambda x: valences[get_element(x)], reverse=True)
                self.ignored_terminal_group_atoms.append(atom)
                self.ignored_terminal_group_atoms.append(sorted_bond_partners[1])
                self.ignored_terminal_group_atoms.append(sorted_bond_partners[2])
                self.ignored_terminal_group_atoms.append(sorted_bond_partners[3])


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
        for i, _ in enumerate(coords2):
            coords2[i] = np.matmul(coords2[i], R)

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
