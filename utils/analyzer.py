import os

import numpy as np
#import multiprocessing as mp

from utils.helper import get_element, valences
from SUPRAConformer.structure import Structure
from typing import Union
from queue import Queue
from scipy.optimize import linear_sum_assignment



class Analyzer:

    def __init__(self):
        pass



    def _count_conformers_dir(self, path: str) -> int:
        path = os.path.abspath(path)
        conformers = os.listdir(path)
        return len(conformers)
    


    def _count_conformers_file(self, path: str) -> int:
        with open(path, "r") as infile:
            n_atoms = int(infile.readline().split()[0])
            conformer_counter = 0
            line_iter = 0
            infile.seek(0)
            for line in infile:
                if (line_iter >= n_atoms+1):
                    line_iter = 0
                    conformer_counter += 1
                line_iter += 1
        return conformer_counter
    


    def check_for_duplicates(
        self, 
        xyz_file: str,
        path: str,
        rmsd_threshold: float=0.1,
        matching: str="normal",
        use_energy: bool = False,
        ignore: str=None,
    ) -> str:
        conformer1 = Structure()
        conformer2 = Structure()
        conformers = os.listdir(path)

        if ignore == "methyl":
            conformer1.get_structure(xyz_file, read_energy=use_energy)
            for atom in self._find_methyl_group_atoms(conformer1.bond_partners):
                del conformer1.coords[atom]
        elif ignore == "terminal":
            conformer1.get_structure(xyz_file, read_energy=use_energy)
            for atom in self._find_terminal_group_atoms(conformer1.bond_partners):
                del conformer1.coords[atom]
        else:
            conformer1.read_xyz(xyz_file, read_energy=use_energy)

        for other_file in conformers:
            other_file = os.path.join(path, other_file)
            if other_file == xyz_file:
                continue
            if ignore == "methyl":
                conformer2.get_structure(other_file, read_energy=use_energy)
                for atom in self._find_methyl_group_atoms(conformer1.bond_partners):
                    del conformer1.coords[atom]
            elif ignore == "terminal":
                conformer2.get_structure(other_file, read_energy=use_energy)
                for atom in self._find_terminal_group_atoms(conformer2.bond_partners):
                    del conformer2.coords[atom]
            else:
                conformer2.read_xyz(other_file, read_energy=use_energy)
            if matching == "loose":
                rmsd = self._calc_rmsd(conformer1.coords, conformer2.coords)
            elif matching == "normal":
                rmsd = self._rmsd(conformer1.coords, conformer2.coords)
            elif matching == "tight":
                if not ignore:
                    conformer1.get_connectivity()
                    conformer2.get_connectivity()
                rmsd = self._rmsd_tight(
                    conformer1.coords,
                    conformer1.bond_partners,
                    conformer2.coords,
                    conformer2.bond_partners
                )
            if rmsd <= rmsd_threshold:
                if (conformer1.energy and conformer2.energy):
                        if conformer1.energy < conformer2.energy:
                            return xyz_file
                        else:
                            return other_file
                else:
                    return xyz_file
        
        return None
                


    def compare_ensembles_dirs(self, path1: str, path2: str, rmsd_threshold: float=0.1, ignore: str=None, matching: str="normal") -> list:
        conformer1 = Structure()
        conformer2 = Structure()

        #path1 = os.path.abspath(path1)
        #path2 = os.path.abspath(path2)
        conformers1 = os.listdir(path1)
        conformers2 = os.listdir(path2)
        if len(conformers1) > len(conformers2):
            conformers1, conformers2 = conformers2, conformers1
            path1, path2 = path2, path1

        counter = 0
        #m = round(len(conformers1)/50)
        #print("[", end="", flush="True")
        for file1 in conformers1:
            #if (i % m) == 0:
            #    print("=", end="", flush=True)
            if ignore == "methyl":
                conformer1.get_structure(os.path.join(path1, file1))
                for atom in self._find_methyl_group_atoms(conformer1.bond_partners):
                    del conformer1.coords[atom]
            elif ignore == "terminal":
                conformer1.get_structure(os.path.join(path1, file1))
                for atom in self._find_terminal_group_atoms(conformer1.bond_partners):
                    del conformer1.coords[atom]
            else:
                conformer1.read_xyz(os.path.join(path1, file1))
            for file2 in conformers2:
                if ignore == "methyl":
                    conformer2.get_structure(os.path.join(path2, file2))
                    for atom in self._find_methyl_group_atoms(conformer2.bond_partners):
                        del conformer2.coords[atom]
                elif ignore == "terminal":
                    conformer2.get_structure(os.path.join(path2, file2))
                    for atom in self._find_terminal_group_atoms(conformer2.bond_partners):
                        del conformer2.coords[atom]
                else:
                    conformer2.read_xyz(os.path.join(path2, file2))
                if matching == "loose":
                    rmsd = self._calc_rmsd(conformer1.coords, conformer2.coords)
                elif matching == "normal":
                    rmsd = self._rmsd(conformer1.coords, conformer2.coords)
                elif matching == "tight":
                    if not ignore:
                        conformer1.get_connectivity()
                        conformer2.get_connectivity()
                    rmsd = self._rmsd_tight(conformer1.coords, conformer1.bond_partners, conformer2.coords, conformer2.bond_partners)
                if rmsd <= rmsd_threshold:
                    counter += 1
                    break
        #print("]")

        return [path1, len(conformers1), path2, len(conformers2), counter]



    def compare_ensembles_files(self, input_path1: str, input_path2: str, rmsd_threshold: float=0.1, ignore: str=None, matching: str="normal") -> list:
        ensemble_file1 = os.path.abspath(input_path1)
        ensemble_file2 = os.path.abspath(input_path2)
        dir_ensemble_file1 = os.path.dirname(ensemble_file1)
        dir_ensemble_file2 = os.path.dirname(ensemble_file2)
        dir_ensemble_file1 = os.path.join(dir_ensemble_file1, "conformers1")
        dir_ensemble_file2 = os.path.join(dir_ensemble_file2, "conformers2")

        os.makedirs(dir_ensemble_file1)
        os.makedirs(dir_ensemble_file2)

        struc_filename = "conformer"

        with open(ensemble_file1, "r") as infile:
            n_atoms = int(infile.readline().split()[0])
            file_counter = 0
            line_iter = 0
            infile.seek(0)
            new_struc = os.path.join(dir_ensemble_file1, f"{struc_filename}{file_counter}.xyz")
            new_struc_file = open(new_struc, "w")
            for line in infile:
                if (line_iter >= n_atoms+2):
                    new_struc_file.close()
                    line_iter = 0
                    file_counter += 1
                    new_struc = os.path.join(dir_ensemble_file1, f"{struc_filename}{file_counter}.xyz")
                    new_struc_file = open(new_struc, "w")
                print(line, end="", file=new_struc_file)
                line_iter += 1
            new_struc_file.close()
        
        with open(ensemble_file2, "r") as infile:
            n_atoms = int(infile.readline().split()[0])
            file_counter = 0
            line_iter = 0
            infile.seek(0)
            new_struc = os.path.join(dir_ensemble_file2, f"{struc_filename}{file_counter}.xyz")
            new_struc_file = open(new_struc, "w")
            for line in infile:
                if (line_iter >= n_atoms+2):
                    new_struc_file.close()
                    line_iter = 0
                    file_counter += 1
                    new_struc = os.path.join(dir_ensemble_file2, f"{struc_filename}{file_counter}.xyz")
                    new_struc_file = open(new_struc, "w")
                print(line, end="", file=new_struc_file)
                line_iter += 1
            new_struc_file.close()
        
        path_less, n_conformers_less, path_more, n_conformers_more, overlap = self.compare_ensemble_dirs(
                path1=dir_ensemble_file1,
                path2=dir_ensemble_file2,
                rmsd_threshold=rmsd_threshold,
                ignore=ignore,
                matching=matching
        )

        if path_more == dir_ensemble_file1 and path_less == dir_ensemble_file2:
            input_path1, input_path2 = input_path2, input_path1
               
        return [input_path1, n_conformers_less, input_path2, n_conformers_more, overlap]




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



    def remove_doubles_dir(self, path: str, rmsd_threshold: float=0.1, ignore: str=None, use_energy: bool = False, matching: str="normal") -> int:
        conformer1 = Structure()
        conformer2 = Structure()
        path = os.path.abspath(path)
        conformers = os.listdir(path)
        #m = len(conformers)/50
        #n = 1
        delete_files = [0 for _ in range(len(conformers))]

        #print(f"{'_'*50}")
        #iprint("|", end="", flush=True)
        for i, file1 in enumerate(conformers):
            #if index > m*n:
            #    print("#", end="", flush=True)
            #    n += 1
            if ignore == "methyl":
                conformer1.get_structure(os.path.join(path, file1), read_energy=use_energy)
                for atom in self._find_methyl_group_atoms(conformer1.bond_partners):
                    del conformer1.coords[atom]
            elif ignore == "terminal":
                conformer1.get_structure(os.path.join(path, file1), read_energy=use_energy)
                for atom in self._find_terminal_group_atoms(conformer1.bond_partners):
                    del conformer1.coords[atom]
            else:
                conformer1.read_xyz(os.path.join(path, file1), read_energy=use_energy)
            for j, file2 in enumerate(conformers[i+1:], start=i+1):
                if ignore == "methyl":
                    conformer2.get_structure(os.path.join(path, file2), read_energy=use_energy)
                    for atom in self._find_methyl_group_atoms(conformer2.bond_partners):
                        del conformer2.coords[atom]
                elif ignore == "terminal": 
                    conformer2.get_structure(os.path.join(path, file2), read_energy=use_energy)
                    for atom in self._find_terminal_group_atoms(conformer2.bond_partners):
                        del conformer2.coords[atom]
                else:
                    conformer2.read_xyz(os.path.join(path, file2), read_energy=use_energy)
                if matching == "loose":
                    rmsd = self._calc_rmsd(conformer1.coords, conformer2.coords)
                elif matching == "normal":
                    rmsd = self._rmsd(conformer1.coords, conformer2.coords)
                elif matching == "tight":
                    if not ignore:
                        conformer1.get_connectivity()
                        conformer2.get_connectivity()
                    rmsd = self._rmsd_tight(conformer1.coords, conformer1.bond_partners, conformer2.coords, conformer2.bond_partners)
                if rmsd <= rmsd_threshold:
                    if (conformer1.energy and conformer2.energy):
                        if conformer1.energy < conformer2.energy:
                            delete_files[j] = 1
                        else:
                            delete_files[i] = 1
                    else:
                        delete_files[j] = 1
        
        #if sum(delete_files):
        #    raw_path = os.path.dirname(path)
        #    out_dir = os.path.join(raw_path, "SUPRA_Output")
        #    os.system(f"mkdir {out_dir}")

        counter = len(conformers)
        #counter = 0
        for i, delete in enumerate(delete_files):
            if delete:
                os.remove(os.path.join(path, conformers[i]))
                counter -= 1
            #if not delete:
            #    os.system(f"cp {os.path.join(path, conformers[i])} {os.path.join(out_dir, conformers[i])}")
            #    counter += 1

        #print("\n")
        
        return counter
    


    def remove_doubles_file(self, ensemble_file: str, rmsd_threshold: float=0.1, ignore: str=None, use_energy: bool = False, matching: str="normal") -> int:
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
        
        n_conformers = self.remove_doubles_dir(workdir, rmsd_threshold, ignore, use_energy, matching)
        
        with open(os.path.join(dir_ensemble_file, "SUPRA_output_ensemble.xyz"), "w") as outfile:
            conformers = os.listdir(workdir)
            for conformer in conformers:
                with open(os.path.join(workdir, conformer), "r") as infile:
                    print(infile.read(), end="", file=outfile)
        
        os.system(f"rm -rf {workdir}")

        return n_conformers



    def _find_methyl_group_atoms(self, structure: dict) -> list:
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



    def _find_terminal_group_atoms(self, structure: dict) -> list:
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
    


    def _rmsd(self, coords1: dict, coords2: dict) -> float:
        if len(coords1.keys()) != len(coords2.keys()):
            return 1000.0
        new_coords1, new_coords2 = self._hungarian_match(coords1, coords2)
        return self._calc_rmsd(new_coords1, new_coords2)
    


    def SCHK(self, mol1: Structure, mol2: Structure) -> tuple:
        if len(mol1.number_of_atoms != mol2.number_of_atoms):
            return 1000.0
        
        # atoms = ...

        elements1 = [get_element(atom) for atom in mol1.coords.keys()]
        elements2 = [get_element(atom) for atom in mol2.coords.keys()]

        new_coords1 = np.array(list(mol1.coords.values()))
        new_coords2 = np.array(list(mol2.coords.values()))

        #for bond in mol1.bonds:

        new_coords1, elements1, new_coords2, elements2 = self._rmsd_hungarian(
            new_coords1, elements1, new_coords2, elements2
        )

        current_rmsd = self._calc_rmsd(new_coords1, new_coords2)
        last_rmsd = 1000000.0
        delta = abs(current_rmsd - last_rmsd)
        threshold = 0.1
        max_iter = 10
        iter = 0

        while (delta > threshold and iter < max_iter):
            print(current_rmsd)
            iter += 1
            last_rmsd = current_rmsd
            kabsch_coords1, kabsch_coords2 = self._kabsch(new_coords1, new_coords2)
            new_coords1, elements1, new_coords2, elements2 = self._rmsd_hungarian(
                kabsch_coords1, elements1, kabsch_coords2, elements2
            )
            current_rmsd = self._calc_rmsd(new_coords1, new_coords2)
            delta = abs(current_rmsd - last_rmsd)
        
        return new_coords1, new_coords2
    


    def _rmsd_tight(self, coords1: dict, connectivity1: dict, coords2: dict, connectivity2: dict,) -> float:
        if len(coords1.keys()) != len(coords2.keys()):
            return 1000.0
        coords1, coords2 = self._match(coords1, connectivity1, coords2, connectivity2)
        kabsch_coords1, kabsch_coords2 = self._kabsch(coords1, coords2)
        return self._calc_rmsd(kabsch_coords1, kabsch_coords2)



    def _kabsch(self, coords1: Union[dict, list, np.array], coords2: Union[dict, list, np.array]) -> tuple:
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
    


    def _calc_rmsd(self, coords1: Union[dict, list, np.array], coords2: Union[dict, list, np.array]) -> float:
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



    def _hungarian_match(self, coords1: dict, coords2: dict) -> tuple:
        elements1 = [get_element(atom) for atom in coords1.keys()]
        elements2 = [get_element(atom) for atom in coords2.keys()]
        n_atoms = len(elements1)
        cost = np.zeros((n_atoms, n_atoms))
        kabsch_coords1, kabsch_coords2 = self._kabsch(coords1, coords2)
        for i in range(n_atoms):
            for j in range(n_atoms):
                if elements1[i] == elements2[j]:
                    element_term = 0.0
                else:
                    element_term = 1000.0
                diff_vec = kabsch_coords1[i] - kabsch_coords2[j]
                cost_value = np.dot(diff_vec, diff_vec) + element_term
                cost[i][j] = cost_value
        row, col = linear_sum_assignment(cost)
        return kabsch_coords1[row], kabsch_coords2[col]
    


    def _rmsd_hungarian(self, coords1: np.array, elements1: list, coords2: np.array, elements2: list) -> list:
        n_atoms = len(elements1)
        cost = np.zeros((n_atoms, n_atoms))
        for i in range(n_atoms):
            for j in range(n_atoms):
                if elements1[i] == elements2[j]:
                    element_term = 0.0
                else:
                    element_term = 1000.0
                diff_vec = coords1[i] - coords2[j]
                cost_value = np.dot(diff_vec, diff_vec) + element_term
                cost[i][j] = cost_value
        row, col = linear_sum_assignment(cost)
        new_elements2 = [elements2[i] for i in col]
        return coords1[row], elements1, coords2[col], new_elements2



    # No. 1 are reference coordinates and structure
    # No. 2 coordinates will be reordered to match No. 1
    def _match(self, coords1: dict, structure1: dict, coords2: dict, structure2: dict) -> tuple:
        # initialization
        pairs = {}
        matched_coords1 = []
        matched_coords2 = []
        list1 = list(coords1.keys())
        list2 = list(coords2.keys())
        spheres1 = {atom: self._spheres(structure1, atom) for atom in list1}
        spheres2 = {atom: self._spheres(structure2, atom) for atom in list2}

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
            elif len(eq) == 1:
                matched_coords1.append(coords1[atom1])
                matched_coords2.append(coords2[eq[0]])
            # security option for cases where no potentially equivalent atom is found,
            # return coordinates with unchanged ordering in this case
            else:
                return list(coords1.values()), list(coords2.values())
        
        # find correct atom pairs from No. 1 and No. 2 via matching of distances to other atoms
        for atom, eq_atoms in pairs.items():
            print(atom)
            d1 = np.zeros(shape=(len(matched_coords1), len(eq_atoms)))
            d2 = np.zeros(shape=(len(matched_coords2), len(eq_atoms)))
            # store distances Pos.(atom No. 1)-Pos.(matched atom i from No. 1) in matrix d1 which is
            # of shape (number of already matched coordinates)x(number of candidate atoms) which means
            # every row stands for a reference atom from No. 1, every column contains the distance of the
            # yet unmatched atom to the respective reference atom
            for i, comp_atom in enumerate(matched_coords1):
                v11 = coords1[atom]
                v12 = comp_atom
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
            print(eq)
            # store coordinates of atom from No. 1 and coordinates of atom from No. 2 with smallest deviation
            matched_coords1.append(coords1[atom])
            matched_coords2.append(coords2[eq])
            # remove candidate from candidate lists
            for candidates in pairs.values():
                if eq in candidates:
                    candidates.remove(eq)

        return np.array(matched_coords1), np.array(matched_coords2)



    def _spheres(self, structure: dict, start: str) -> list:
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
    