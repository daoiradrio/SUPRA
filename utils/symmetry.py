import numpy as np

from SUPRAConformer.structure import Structure
from utils.rotationaxis import RotationAxis
from utils.helper import atom_in_torsions, get_element



class Symmetry:

    def __init__(self):
        self.possible_orders = [
            #6, # corresponds to angle increment  60°
            #4, # corresponds to angle increment  90°
            #3, # corresponds to angle increment 120°
            2  # corresponds to angle increment 180°
        ]


    def _get_torsion_group(self, torsions: list, connectivity: dict, atom: str, last_atom: str, status: dict, torsion_atoms: list):
        if status[atom] == "SEEN":
            return
        else:
            status[atom] = "SEEN"
            for neighbor in connectivity[atom]:
                if neighbor != last_atom:
                    if not atom_in_torsions(torsions, neighbor):
                        torsion_atoms.append(neighbor)
                        self._get_torsion_group(torsions, connectivity, neighbor, atom, status, torsion_atoms)
        return


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
            print(f"Side of {atom1}: {left_rot_sym}")
            print(f"Side of {atom2}: {right_rot_sym}")
            print()

    

    def rot_sym_along_bond(self, mol: Structure, rot_axis: RotationAxis, rot_atoms: list, order: int) -> bool:
        atom_used = [0 for i in range(mol.number_of_atoms)]

        for i, atom_i in enumerate(rot_atoms):
            if (atom_used[i]):
                continue
            coords_i = mol.coords[atom_i]
            best_j = i
            best_distance = 1.0
            symmetric_coords = rot_axis.rotate_atom(coords_i, 360.0/float(order))
            for j, atom_j in enumerate(rot_atoms[i+1:], start=i+1):
                if (atom_used[j]):
                    continue
                if (get_element(atom_i) != get_element(atom_j)):
                    continue
                coords_j = mol.coords[atom_j]
                distance = np.linalg.norm(coords_j - symmetric_coords)
                if (distance < best_distance):
                    best_j = j
                    best_distance = distance
            if (best_distance > 0.5):
                return False
            atom_used[best_j] = 1
        
        return True
    

    def rot_order_along_bond(self, mol: Structure, rot_atoms: list, from_coords: np.array, to_coords: np.array) -> int:
        rot_axis = RotationAxis(from_coords, to_coords)
        for order in self.possible_orders:
            if (self.rot_sym_along_bond(mol, rot_axis, rot_atoms, order)):
            #if (self.rot_sym_along_bond(mol, from_coords, to_coords, rot_atoms, order)):
                return order
        return 1