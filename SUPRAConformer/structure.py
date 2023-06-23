import os # path.exists

import numpy as np # sqrt

from utils.bond import Bond
from utils.helper import covalence_radii_single, covalence_radii_double, covalence_radii_triple, get_element



class Structure:

    def __init__(self, file: str = None):
        self.number_of_atoms = 0
        self.coords = {}
        self.bond_partners = {}
        self.bond_orders = {}
        self.bonds = []
        self.energy = 0
        if file:
            self.get_structure(file)



    def get_structure(self, filename: str, read_energy: bool = False):
        if not filename:
            filename = input("Path of the .xyz-File: ")
        if not os.path.exists(filename):
            print("STRUCTURE MODULE: File not found at given path.")
            return
        self.read_xyz(filename, read_energy)
        self.get_connectivity()



    def read_xyz(self, filename: str, read_energy: bool = False):
        self.coords = {}
        with open(filename, "r") as input_file:
            for i, line in enumerate(input_file):
                if i == 0:
                    self.number_of_atoms = int(line)
                elif read_energy and i == 1:
                    self.energy = float(line.split()[-1])
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
            #max_valence = valences[self.get_element(atom1)]
            valence = 0
            for atom2 in atoms[i+1:]:
                #if valence == max_valence:
                #    break
                coords2 = self.coords[atom2]
                bond_order = self._check_connectivity(atom1, coords1, atom2, coords2)
                if bond_order:
                    new_bond = Bond()
                    new_bond.atom1 = atom1
                    new_bond.atom2 = atom2
                    new_bond.bond_order = bond_order
                    self.bonds.append(new_bond)
                    #self.bonds.append((atom1, atom2))
                    #self.bonds.append(sorted([atom1, atom2], key=lambda label: get_number(label)))
                    self.bond_partners[atom1].append(atom2)
                    self.bond_partners[atom2].append(atom1)
                    #self.bond_orders[(atom1, atom2)] = bond_order
                    valence += 1



    # calculate connectivity and bond order
    def _check_connectivity(self, atom1: str, coord1: np.array, atom2: str, coord2: np.array) -> int:
        # initialize variables
        tolerance = 0.08
        bond_order = 0
        element1 = get_element(atom1)
        element2 = get_element(atom2)
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
