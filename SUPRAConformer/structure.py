# covalence radii, max. valences, atom label to element symbol converter
from utils.helper import valences, covalence_radii_single, covalence_radii_double, covalence_radii_triple
import numpy as np # sqrt
import os # path.exists



class Structure:

    def __init__(self, file: str = None):
        self.number_of_atoms = 0
        self.coords = {}
        self.bond_partners = {}
        self.bond_orders = {}
        self.bonds = []
        if file:
            self.get_structure(file)


    # translation of atom label into element symbol
    @staticmethod
    def get_element(label: str) -> str:
        # if letter on second index: two letter element symbol, return first two chars of label
        if label[1].isalpha():
            return (label[0] + label[1])
        # if no letter on second index: one letter element symbol, return first char of label
        else:
            return label[0]

    
    @staticmethod
    def get_number(label: str) -> int:
        if label[-2].isdigit():
            return int(label[-2] + label[-1])
        else:
            return int(label[-1])


    def get_structure(self, filename: str):
        if not filename:
            filename = input("Path of the .xyz-File: ")
        if not os.path.exists(filename):
            print("STRUCTURE MODULE: File not found at given path.")
            return
        self.read_xyz(filename)
        self.get_connectivity()


    def read_xyz(self, filename: str):
        self.coords = {}
        with open(filename, "r") as input_file:
            for i, line in enumerate(input_file):
                if i == 0:
                    self.number_of_atoms = int(line)
                elif i >= 2:
                    element, x, y, z = line.split()
                    self.coords[f"{element}{i-2}"] = np.array([float(x), float(y), float(z)])


    def get_connectivity(self):
        atoms = list(self.coords.keys())
        self.bond_partners = {atom: [] for atom in atoms}
        self.bond_orders = {}
        self.bonds = []
        for i, atom1 in enumerate(atoms):
            coords1 = self.coords[atom1]
            max_valence = valences[self.get_element(atom1)]
            valence = 0
            for atom2 in atoms[i+1:]:
                if valence == max_valence:
                    break
                else:
                    coords2 = self.coords[atom2]
                    bond_order = self._check_connectivity(atom1, coords1, atom2, coords2)
                    if bond_order:
                        self.bonds.append((atom1, atom2))
                        self.bond_partners[atom1].append(atom2)
                        self.bond_partners[atom2].append(atom1)
                        self.bond_orders[(atom1, atom2)] = bond_order
                        valence += 1


    # TODO: ÜBERPRÜFUNG HINZUFÜGEN (PFAD KORREKT, DATEI LESBAR, ETC.)
    # read atom coordinates from .xyz-file, calculate connectivity
    def old_get_structure(self, filename: str) -> None:
        if not filename:
            filename = input("Path of the .xyz-File: ")
        if not os.path.exists(filename):
            print("STRUCTURE MODULE: File not found at given path.")
            return
        self.coords = {}
        self.bond_partners = {}
        self.bond_orders = {}
        self.bonds = []
        with open(filename, "r") as input_file:
            for i, line in enumerate(input_file):
                if i == 0:
                    self.number_of_atoms = int(line)
                elif i >= 2:
                    element, x, y, z = line.split()
                    self.coords[f"{element}{i}"] = np.array([float(x), float(y), float(z)])
                    self.bond_partners[f"{element}{i}"] = []
        atoms = list(self.coords.keys())
        for i, atom1 in enumerate(atoms):
            coords1 = self.coords[atom1]
            max_valence = valences[self.get_element(atom1)]
            valence = 0
            for atom2 in atoms[i+1:]:
                if valence == max_valence:
                    break
                else:
                    coords2 = self.coords[atom2]
                    bond_order = self._check_connectivity(atom1, coords1, atom2, coords2)
                    if bond_order:
                        self.bonds.append((atom1, atom2))
                        self.bond_partners[atom1].append(atom2)
                        self.bond_partners[atom2].append(atom1)
                        self.bond_orders[(atom1, atom2)] = bond_order
                        valence += 1


    # calculate connectivity and bond order
    def _check_connectivity(self, atom1: str, coord1: np.array, atom2: str, coord2: np.array) -> int:
        # initialize variables
        tolerance = 0.08
        bond_order = 0
        element1 = self.get_element(atom1)
        element2 = self.get_element(atom2)
        # distance of atom1 and atom2
        distance = np.linalg.norm(coord1 - coord2)
        # in the following -1000 is returned if elements are not implemented yet
        # single bond length for element of atom1 with element of atom2
        single_bond = covalence_radii_single.get(element1, -1000) + covalence_radii_single.get(element2, -1000)
        # double bond length for element of atom1 with element of atom2
        double_bond = covalence_radii_double.get(element1, -1000) + covalence_radii_double.get(element2, -1000)
        # triple bond length for element of atom1 with element of atom2
        triple_bond = covalence_radii_triple.get(element1, -1000) + covalence_radii_triple.get(element2, -1000)
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
