import os

import numpy as np
import multiprocessing as mp

from utils.helper import get_element, valences
from SUPRAConformer.structure import Structure
from typing import Union
from queue import Queue
from scipy.optimize import linear_sum_assignment



class Analyzer:

    def __init__(self):
        pass



    def compare_structure_sets(self, path1: str, path2: str, rmsd_threshold: float=0.1, ignore: str=None) -> int:
        conformer1 = Structure()
        conformer2 = Structure()

        path1 = os.path.abspath(path1)
        path2 = os.path.abspath(path2)
        conformers1 = os.listdir(path1)
        conformers2 = os.listdir(path2)
        if len(conformers1) > len(conformers2):
            conformers1, conformers2 = conformers2, conformers1
            path1, path2 = path2, path1

        counter = 0
        #m = round(len(conformers1)/50)
        print()
        print("Comparing structures...")
        #print("[", end="", flush="True")
        for i, file1 in enumerate(conformers1):
            #if (i % m) == 0:
            #    print("=", end="", flush=True)
            if ignore == "methyl":
                conformer1.get_structure(os.path.join(path1, file1))
                for atom in self.get_methyl_group_atoms(conformer1.bond_partners):
                    print(atom)
                    del conformer1.coords[atom]
            elif ignore == "all":
                conformer1.get_structure(os.path.join(path1, file1))
                for atom in self.get_terminal_group_atoms(conformer1.bond_partners):
                    del conformer1.coords[atom]
            else:
                conformer1.read_xyz(os.path.join(path1, file1))
            for file2 in conformers2:
                if ignore == "methyl":
                    conformer2.get_structure(os.path.join(path2, file2))
                    for atom in self.get_methyl_group_atoms(conformer2.bond_partners):
                        del conformer2.coords[atom]
                elif ignore == "all":
                    conformer2.get_structure(os.path.join(path2, file2))
                    for atom in self.get_terminal_group_atoms(conformer2.bond_partners):
                        del conformer2.coords[atom]
                else:
                    conformer2.read_xyz(os.path.join(path2, file2))
                if self.rmsd(conformer1.coords, conformer2.coords) <= rmsd_threshold:
                    counter += 1
                    break
        #print("]")
        print("Comparing structures done.")
        print()
        print(f"Path 1: {path1}")
        print(f"Path 2: {path2}")
        print()
        print(f"Number of structures in Path 1: {len(conformers1)}")
        print(f"Number of structures in Path 2: {len(conformers2)}")
        print(f"Number of structures of Path 1 in Path 2: {counter}")
        print()

        return counter



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



    def remove_doubles(self, path: str, rmsd_threshold: float=0.1, ignore: str=None, use_energy: bool = False) -> int:
        conformer1 = Structure()
        conformer2 = Structure()
        path = os.path.abspath(path)
        conformers = os.listdir(path)
        #m = len(conformers)/50
        #n = 1
        counter = len(conformers)
        delete_files = [0 for i in range(counter)]

        #print(f"{'_'*50}")
        #iprint("|", end="", flush=True)
        for i, file1 in enumerate(conformers):
            #if index > m*n:
            #    print("#", end="", flush=True)
            #    n += 1
            if ignore == "methyl":
                conformer1.get_structure(os.path.join(path, file1), read_energy=use_energy)
                for atom in self.get_methyl_group_atoms(conformer1.bond_partners):
                    del conformer1.coords[atom]
            elif ignore == "all":
                conformer1.get_structure(os.path.join(path, file1), read_energy=use_energy)
                for atom in self.get_terminal_group_atoms(conformer1.bond_partners):
                    del conformer1.coords[atom]
            else:
                conformer1.read_xyz(os.path.join(path, file1), read_energy=use_energy)
            for j, file2 in enumerate(conformers[i+1:], start=i+1):
                if ignore == "methyl":
                    conformer2.get_structure(os.path.join(path, file2), read_energy=use_energy)
                    for atom in self.get_methyl_group_atoms(conformer2.bond_partners):
                        del conformer2.coords[atom]
                elif ignore == "all": 
                    conformer2.get_structure(os.path.join(path, file2), read_energy=use_energy)
                    for atom in self.get_terminal_group_atoms(conformer2.bond_partners):
                        del conformer2.coords[atom]
                else:
                    conformer2.read_xyz(os.path.join(path, file2), read_energy=use_energy)
                if self.rmsd(conformer1.coords, conformer2.coords) <= rmsd_threshold:
                    if (conformer1.energy and conformer2.energy):
                        if conformer1.energy < conformer2.energy:
                            delete_files[j] = 1
                        else:
                            delete_files[i] = 1
                    else:
                        delete_files[j] = 1

        for i, delete in enumerate(delete_files):
            if delete:
                os.remove(os.path.join(path, conformers[i]))
                counter -= 1

        #print("\n")
        
        return counter
    


    def remove_doubles_ensemble_file(self, ensemble_file: str, rmsd_threshold: float=0.1, ignore: str=None, use_energy: bool = False):
        ensemble_file = os.path.abspath(ensemble_file)
        dir_ensemble_file = os.path.dirname(ensemble_file)
        workdir = os.path.join(dir_ensemble_file, "conformers")
        struc_filename = "conformer"

        os.makedirs(workdir)

        with open(ensemble_file, "r") as infile:
            n_atoms = int(infile.readline().split()[0])
            file_counter = 0
            line_iter = 0
            infile.seek(0)
            new_struc = os.path.join(workdir, f"{struc_filename}{file_counter}.xyz")
            new_struc_file = open(new_struc, "w")
            for line in infile:
                if (line_iter >= n_atoms+2):
                    new_struc_file.close()
                    line_iter = 0
                    file_counter += 1
                    new_struc = os.path.join(workdir, f"{struc_filename}{file_counter}.xyz")
                    new_struc_file = open(new_struc, "w")
                print(line, end="", file=new_struc_file)
                line_iter += 1
            new_struc_file.close()
        
        self.remove_doubles(workdir, rmsd_threshold, ignore, use_energy)
        
        with open(os.path.join(dir_ensemble_file, "supra_ensemble.xyz"), "w") as outfile:
            conformers = os.listdir(workdir)
            for conformer in conformers:
                with open(os.path.join(workdir, conformer), "r") as infile:
                    print(infile.read(), end="", file=outfile)
        
        os.system(f"rm -rf {workdir}")



    def old_remove_doubles(self, path: str, rmsd_threshold: float=0.1, ignore: str=None) -> int:
        print("Performing removal of duplicate structures...")

        conformer1 = Structure()
        conformer2 = Structure()
        path = os.path.abspath(path)
        conformers = os.listdir(path)
        #m = len(conformers)/50
        #n = 1
        counter = len(conformers)

        #print(f"{'_'*50}")
        #iprint("|", end="", flush=True)
        for index, file1 in enumerate(conformers):
            #if index > m*n:
            #    print("#", end="", flush=True)
            #    n += 1
            if ignore == "methyl":
                conformer1.get_structure(os.path.join(path, file1))
                for atom in self.get_methyl_group_atoms(conformer1.bond_partners):
                    del conformer1.coords[atom]
            elif ignore == "all":
                conformer1.get_structure(os.path.join(path, file1))
                for atom in self.get_terminal_group_atoms(conformer1.bond_partners):
                    del conformer1.coords[atom]
            else:
                conformer1.read_xyz(os.path.join(path, file1))
            for file2 in conformers[index + 1:]:
                if ignore == "methyl":
                    conformer2.get_structure(os.path.join(path, file2))
                    for atom in self.get_methyl_group_atoms(conformer2.bond_partners):
                        del conformer2.coords[atom]
                elif ignore == "all": 
                    conformer2.get_structure(os.path.join(path, file2))
                    for atom in self.get_terminal_group_atoms(conformer2.bond_partners):
                        del conformer2.coords[atom]
                else:
                    conformer2.read_xyz(os.path.join(path, file2))
                if self.rmsd(conformer1.coords, conformer2.coords, rmsd_threshold) <= rmsd_threshold:
                    os.remove(os.path.join(path, file1))
                    counter -= 1
                    break
        #print("\n")
        print("Removal of double structures done.")
        print(f"Individual conformers in {path}: {counter}")
   
        return counter
    


    def get_methyl_group_atoms(self, structure: dict) -> list:
        methyl_group_atoms = []
        for atom, bond_partners in structure.items():
            if get_element(atom) == "C":
                terminal_count = 0
                bond_partners = [get_element(bond_partner) for bond_partner in bond_partners]
                terminal_count += bond_partners.count("H")
                terminal_count += bond_partners.count("F")
                terminal_count += bond_partners.count("Cl")
                terminal_count += bond_partners.count("Br")
                terminal_count += bond_partners.count("I")
                if terminal_count == valences[get_element(atom)]-1:
                    sorted_bond_partners = sorted(
                                            structure[atom],
                                            key=lambda x: valences[get_element(x)],
                                            reverse=True
                                           )
                    methyl_group_atoms.append(atom)
                    methyl_group_atoms.append(sorted_bond_partners[1])
                    methyl_group_atoms.append(sorted_bond_partners[2])
                    methyl_group_atoms.append(sorted_bond_partners[3])
        return methyl_group_atoms



    def get_terminal_group_atoms(self, structure: dict) -> list:
        terminal_group_atoms = []
        for atom, bond_partners in structure.items():
            if not get_element(atom) in ["H", "F", "Cl", "Br", "I"]:
                terminal_count = 0
                bond_partners = [get_element(bond_partner) for bond_partner in bond_partners]
                terminal_count += bond_partners.count("H")
                terminal_count += bond_partners.count("F")
                terminal_count += bond_partners.count("Cl")
                terminal_count += bond_partners.count("Br")
                terminal_count += bond_partners.count("I")
                if terminal_count == valences[get_element(atom)]-1:
                    sorted_bond_partners = sorted(
                                            structure[atom],
                                            key=lambda x: valences[get_element(x)],
                                            reverse=True
                                           )
                    terminal_group_atoms.append(atom)
                    for terminal_atom in sorted_bond_partners[1:]:
                        terminal_group_atoms.append(terminal_atom)
        return terminal_group_atoms
    


    def rmsd(self, coords1: dict, coords2: dict) -> float:
        if len(coords1.keys()) != len(coords2.keys()):
            return 1000.0
        elements1 = [get_element(atom) for atom in coords1.keys()]
        elements2 = [get_element(atom) for atom in coords2.keys()]
        n_atoms = len(elements1)
        cost = np.zeros((n_atoms, n_atoms))
        kabsch_coords1, kabsch_coords2 = self.kabsch(coords1, coords2)
        for i in range(n_atoms):
            for j in range(i+1):
                if elements1[i] == elements2[j]:
                    element_term = 0.0
                else:
                    element_term = 1000.0
                diff_vec = kabsch_coords1[i] - kabsch_coords2[j]
                cost_value = np.dot(diff_vec, diff_vec) + element_term
                cost[i][j] = cost_value
                cost[j][i] = cost_value
        row, col = linear_sum_assignment(cost)
        return self.calc_rmsd(kabsch_coords1[row], kabsch_coords2[col])
    


    def rmsd_tight(self, coords1: dict, connectivity1: dict, coords2: dict, connectivity2: dict,) -> float:
        coords1, coords2 = self.match(coords1, connectivity1, coords2, connectivity2)
        return self.calc_rmsd(coords1, coords2)



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
    


    def calc_rmsd(self, coords1: Union[dict, list, np.array], coords2: Union[dict, list, np.array]) -> float:
        if type(coords1) == dict:
            coords1 = list(coords1.values())
        if type(coords2) is dict:
            coords2 = list(coords2.values())
        n = len(coords1)
        delta_sum = 0.0
        for i in range(n):
            delta_sum += ((coords1[i][0] - coords2[i][0])**2 + \
                          (coords1[i][1] - coords2[i][1])**2 + \
                          (coords1[i][2] - coords2[i][2])**2)
        return np.sqrt((1.0/float(n)) * delta_sum)



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



