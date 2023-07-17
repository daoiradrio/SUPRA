import numpy as np

from SUPRAConformer.structure import Structure
from utils.rotationaxis import RotationAxis
from utils.helper import get_element



class Symmetry:

    def __init__(self):
        self.possible_rot_orders = [
            6, # corresponds to angle increment  60°
            4, # corresponds to angle increment  90°
            3, # corresponds to angle increment 120°
            2  # corresponds to angle increment 180°
        ]



    def check_rot_sym_of_torsions(self, mol: Structure, torsions: list) -> None:
        for atom1, atom2 in torsions:
            # rotational symmetry left side of torsion bond
            status = {atom: "UNKNOWN" for atom in mol.coords.keys()}
            left_torsion_atoms = []
            self._get_torsion_group(torsions, mol.bond_partners, atom1, atom2, status, left_torsion_atoms)
            left_rot_sym = self.rot_order_along_bond(mol, left_torsion_atoms, mol.coords[atom1], mol.coords[atom2])
            # rotational symmetry right side of torsion bond
            status = {atom: "UNKNOWN" for atom in mol.coords.keys()}
            right_torsion_atoms = []
            self._get_torsion_group(torsions, mol.bond_partners, atom2, atom1, status, right_torsion_atoms)
            right_rot_sym = self.rot_order_along_bond(mol, right_torsion_atoms, mol.coords[atom1], mol.coords[atom2])

            print(f"{atom1} {atom2}")
            print(f"{left_torsion_atoms}")
            print(f"Side of {atom1}: {left_rot_sym}")
            print(f"{right_torsion_atoms}")
            print(f"Side of {atom2}: {right_rot_sym}")
            print()

    

    def rot_sym_along_bond(self, mol: Structure, rot_axis: RotationAxis, rot_atoms: list, order: int) -> bool:
        for atom_i in rot_atoms:
            coords_i = mol.coords[atom_i]
            best_distance = 1.0
            symmetric_coords = rot_axis.rotate(coords_i, 360.0/float(order))
            for atom_j in rot_atoms:
                if (atom_i == atom_j):
                    continue
                if (get_element(atom_i) != get_element(atom_j)):
                    continue
                coords_j = mol.coords[atom_j]
                distance = np.linalg.norm(coords_j - symmetric_coords)
                if (distance < best_distance):
                    best_distance = distance
            if (best_distance > 0.5):
                return False
        return True
    

    def rot_order_along_bond(self, mol: Structure, rot_atoms: list, from_coords: np.array, to_coords: np.array) -> int:
        rot_axis = RotationAxis(from_coords, to_coords)
        for order in self.possible_orders:
            if (self.rot_sym_along_bond(mol, rot_axis, rot_atoms, order)):
                return order
        return 1