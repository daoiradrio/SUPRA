from Helper import covalence_radii, valences, get_element # covalence radii, max. valences,
                                                          # atom label to element symbol converter
import numpy as np # sqrt
import pandas as pd # reading atom coordinates from .xyz-file (csv-file)
import os # path.exists
import queue



class Structure:

    def __init__(self, file: str=None):
        # coordinates of all atoms in input structure
        # coords = {E0:[x0, y0, z0], E1:[x1, y1, z1], ... }
        self.coords = dict()

        # connectivity of all atoms in input structure (undirected graph)
        # structure = {E0:[E1, E2, E3, ... ], E1:[E0, E4, ], ... }
        self.structure = dict()

        # nonterminal single bonds, unfiltered rotatable bonds
        # torsions = [[E0, E1], [E2, E3], ... ]
        self.torsions = list()

        self.torsion_atoms = list()

        if file:
            #HIER ENTSTEHEN PROBLEME WENN ORCA-KOORDINATEN REINGEGEBEN WERDEN (ZWEI ZUSÄTZLICHE KOPFZEILEN)
            self.get_structure(file)


    # WOFÜR ÜBERHAUPT CHECK AUF DREIFACHBINDUNG??
    #
    # calculate connectivity and bond order
    def check_connectivity(self, atom1: str, coord1: list, atom2: str, coord2: list, terminal: bool) -> int:
        # initialize variables
        tolerance = 0.08
        bond_order = 0

        # distance of atom1 and atom2
        distance = np.sqrt((coord1[0] - coord2[0])**2 + (coord1[1] - coord2[1])**2 + (coord1[2] - coord2[2])**2)

        # single bond length for element of atom1 with element of atom2
        single_bond = covalence_radii[atom1][0] + covalence_radii[atom2][0]

        # if atom1 or atom2 is terminal atom only checking for single bond necessary
        if terminal:
            if distance <= (single_bond + tolerance):
                bond_order = 1

        # if atom1 and atom2 are chain atoms also checking for multiple bonds
        else:
            # double bond length for element of atom1 with element of atom2
            double_bond = covalence_radii.get(atom1)[1] + covalence_radii.get(atom2)[1]

            # triple bond length for element of atom1 with element of atom2
            triple_bond = covalence_radii.get(atom1)[2] + covalence_radii.get(atom2)[2]

            # checking for triple bond
            if distance <= (triple_bond + tolerance):
                bond_order = 3

            # checking for double bond
            elif distance <= (double_bond + tolerance):
                bond_order = 2

            # checking for single bond
            elif distance <= (single_bond + tolerance):
                bond_order = 1

        # bond order 1 or higher if bond between atom1 and atom2 exist, otherwise 0
        return bond_order


    # read atom coordinates from .xyz-file, calculate connectivity
    def get_structure(self, filename: str):
        if not filename:
            filename = input("Path of the .xyz-File: ")
        if not os.path.exists(filename):
            print("STRUCTURE MODULE: File not found at given path.")
            return -1

        # clean up coordinate and connectivity storage in case of old data
        if self.coords:
            self.coords = dict()
        if self.structure:
            self.structure = dict()

        # ÜBERPRÜFUNG HINZUFÜGEN (PFAD KORREKT, DATEI LESBAR, ETC.)
        # read data from file at given path
        file = pd.read_csv(filename, header=None, delim_whitespace=True)

        # convert data to list
        list = file.values.tolist()

        # loop over every atom, count them
        for counter, line in enumerate(list):
            # get element (first column)
            element = line[0]
            # get coordinates (all other columns)
            coord = line[1:]

            # create label as key for atom (element + number) in coords-dictionary, store coordinates
            # of atom at this key
            self.coords[element + str(counter)] = coord

            # create label as key for atom (element + number) in structure-dictionary
            self.structure[element + str(counter)] = []

        # loop over every key/label in structure-dictionary
        for atom1 in self.structure:
            # get max. number of valence partners for element of current key
            max_valence = valences[get_element(atom1)]

            # for every key in structure-dictionary loop over every key/label in coords-dictionary
            for atom2 in self.coords:
                # initialize flag for terminal atoms
                terminal = False

                # update valence of current atom
                valence = len(self.structure[atom1])

                # if max. valence of current atom is reached jump to next cycle of outer loop
                if valence == max_valence:
                    break

                # avoid checking connectivity of atom with itself, jump to next cycle of inner loop
                elif atom1 == atom2:
                    continue

                else:
                    # if current atom of outer or inner loop is terminal raise flag
                    if get_element(atom1) in ["H", "F", "Cl", "Br", "I"] or \
                       get_element(atom2) in ["H", "F", "Cl", "Br", "I"]:
                        terminal = True

                    # calculate bond order/connectivity
                    bond_order = self.check_connectivity(get_element(atom1), self.coords[atom1],
                                                         get_element(atom2), self.coords[atom2],
                                                         terminal)

                    # if bond between atom of outer loop and atom of inner loop exists store label of inner loop at
                    # key/label of atom of outer loop in structure-dictionary
                    if bond_order:
                        self.structure[atom1].append(atom2)

                        # if bond is nonterminal (flag not raised) single bond store it in list of rotatable bonds
                        if not terminal and bond_order == 1:
                            # store every rotatable bond only once
                            if not atom1 in self.structure[atom2]:
                                self.torsions.append([atom1, atom2])


    # read atom coordinates from .xyz-file written by ORCA (atom number in first line, text in second line),
    # calculate connectivity
    def get_structure_ORCA(self, filename: str):
        if not filename:
            filename = input("Path of the .xyz-File: ")
            if not os.path.exists(filename):
                print("File not found at given path.")
                return -1

        # clean up coordinate and connectivity storage in case of old data
        if self.coords:
            self.coords = dict()
        if self.structure:
            self.structure = dict()

        # read data row wise from file at given path
        with open(filename, "r") as file:
            data = file.read().splitlines()

        # loop over every atom, count them
        for counter, line in enumerate(data[2:]):
            # create label for atom (element + number)
            key = line.split()[0] + str(counter)

            # declare label as key and store atom coordinates in coords-dictionary at key
            self.coords[key] = [float(line.split()[1]), float(line.split()[2]), float(line.split()[3])]

            # declare label as key in structure-dictionary
            self.structure[key] = []

        # loop over every key/label in structure-dictionary
        for atom1 in self.structure:
            # get max. number of valence partners for element of current key
            max_valence = valences[get_element(atom1)]

            # for every key/label in structure-dictionary loop over every key/label in coords-dictionary
            for atom2 in self.coords:
                # initialize flag for terminal atoms
                terminal = False

                # update valence of current atom
                valence = len(self.structure[atom1])

                # if max. valence of current atom is reached jump to next cycle of outer loop
                if valence == max_valence:
                    break

                # avoid checking connectivity of atom with itself, jump to next cycle of inner loop
                elif atom1 == atom2:
                    continue

                else:
                    # if current atom of outer or inner loop is terminal raise flag
                    if get_element(atom1) in ["H", "F", "Cl", "Br", "I"] or \
                       get_element(atom2) in ["H", "F", "Cl", "Br", "I"]:
                        terminal = True

                    # calculate bond order/connectivity
                    bond_order = self.check_connectivity(get_element(atom1), self.coords[atom1],
                                                         get_element(atom2), self.coords[atom2],
                                                         terminal)

                    # if bond between atom of outer loop and atom of inner loop exists store label of inner
                    # loop at key/label of atom of outer loop in structure-dictionary
                    if bond_order:
                        self.structure[atom1].append(atom2)

                        # if bond is nonterminal (flag not raised) single bond store it in list of rotatable bonds
                        if not terminal and bond_order == 1:
                            # store every rotatable bond only once
                            if not atom1 in self.structure[atom2]:
                                self.torsions.append([atom1, atom2])


    def _clashes(self, coords: dict=None):
        if not coords:
            coords = self.coords.copy()

        atoms_list = list()
        for atom in coords:
            atoms_list.append(atom)

        for index, atom1 in enumerate(atoms_list):
            for atom2 in atoms_list[index + 1:]:
                if atom1 != atom2:

                    distance = np.sqrt((coords[atom1][0] - coords[atom2][0]) ** 2 + \
                                       (coords[atom1][1] - coords[atom2][1]) ** 2 + \
                                       (coords[atom1][2] - coords[atom2][2]) ** 2)

                    minimum_distance = covalence_radii[get_element(atom1)][0] + \
                                       covalence_radii[get_element(atom2)][0] + 0.08

                    if distance < minimum_distance:
                        if self._distant(atom1, atom2):
                            return True
        return False


    def _distant(self, atom1: str, atom2: str):
        neighbors = queue.Queue()
        distance = 1
        for neighbor in self.structure[atom1]:
            neighbors.put(neighbor)
        while distance < 4 and neighbors:
            number_neighbors = neighbors.qsize()
            while number_neighbors:
                curr = neighbors.get()
                if curr == atom2:
                    return False
                else:
                    for neighbor in self.structure[curr]:
                        neighbors.put(neighbor)
                number_neighbors -= 1
            distance += 1
        return True
