import numpy as np

from SUPRAConformer.structure import Structure
from utils.rotationaxis import RotationAxis
from utils.helper import get_element, get_number, atom_in_torsions



class Symmetry:

    def __init__(self):
        self._possible_rot_orders = [
            6, # corresponds to angle increment  60°
            4, # corresponds to angle increment  90°
            3, # corresponds to angle increment 120°
            2  # corresponds to angle increment 180°
        ]



    def _find_rot_sym_of_torsions(self, mol: Structure, torsions: list, angle_increment: int) -> None:
        torsion_done = [0 for _ in torsions]

        for i, torsion in enumerate(torsions):
            atom1 = torsion.atom1
            atom2 = torsion.atom2
            torsion.rot_angles = [0]
            # rotational symmetry left side of torsion bond
            status = {atom: "UNKNOWN" for atom in mol.coords.keys()}
            self._find_torsion_group(mol.bond_partners, torsions, atom1, atom2, status, torsion.sym_rot_atoms1)
            torsion.sym_rot_atoms1 = sorted(torsion.sym_rot_atoms1, key=lambda label: get_number(label))
            torsion.rot_sym1 = self.rot_order_along_bond(mol, torsion.sym_rot_atoms1, mol.coords[atom1], mol.coords[atom2])
            # rotational symmetry right side of torsion bond
            status = {atom: "UNKNOWN" for atom in mol.coords.keys()}
            self._find_torsion_group(mol.bond_partners, torsions, atom2, atom1, status, torsion.sym_rot_atoms2)
            torsion.sym_rot_atoms2 = sorted(torsion.sym_rot_atoms2, key=lambda label: get_number(label))
            torsion.rot_sym2 = self.rot_order_along_bond(mol, torsion.sym_rot_atoms2, mol.coords[atom1], mol.coords[atom2])
            # assign rotation angles if already possible
            if (torsion.rot_sym1 > 1)                                                                         and \
               (max(360/torsion.rot_sym1, angle_increment) % min(360/torsion.rot_sym1, angle_increment) == 0) and \
               (torsion.sym_rot_atoms1 == torsion.rot_atoms1):
                # IN DER PRAXIS (WENN NICHT INKREMENT 60 VERWENDET WIRD) KANN SICH DIESE ABFRAGE VERMUTLICH GESPART WERDEN
                if (max(360/torsion.rot_sym2, angle_increment) % min(360/torsion.rot_sym2, angle_increment) == 0) and \
                   (torsion.sym_rot_atoms2 == torsion.rot_atoms2)                                                 and \
                   (torsion.rot_sym1 < torsion.rot_sym2):
                    # case 3
                    j = 1
                    while (j*angle_increment < 360/torsion.rot_sym2):
                        torsion.rot_angles.append(j*angle_increment)
                        j += 1
                    torsion_done[i] = 1
                    continue
                else:
                    # case 2
                    j = 1
                    while (j*angle_increment < 360/torsion.rot_sym1):
                        torsion.rot_angles.append(j*angle_increment)
                        j += 1
                    torsion_done[i] = 1
                    continue
            elif (torsion.rot_sym2 > 1)                                                                         and \
                 (max(360/torsion.rot_sym2, angle_increment) % min(360/torsion.rot_sym2, angle_increment) == 0) and \
                 (torsion.sym_rot_atoms2 == torsion.rot_atoms2):
                # case 3
                j = 1
                while (j*angle_increment < 360/torsion.rot_sym2):
                    torsion.rot_angles.append(j*angle_increment)
                    j += 1
                torsion_done[i] = 1
                continue
            elif (torsion.rot_sym1 == 1 or max(360/torsion.rot_sym1, angle_increment) % min(360/torsion.rot_sym1, angle_increment) != 0) and \
                 (torsion.rot_sym2 == 1 or max(360/torsion.rot_sym2, angle_increment) % min(360/torsion.rot_sym2, angle_increment) != 0):
                # case 1
                for j in range(1, 360//angle_increment):
                    torsion.rot_angles.append(j*angle_increment)
                torsion_done[i] = 1
                continue
            """
            if (torsion.rot_sym1 == 1 or max(360/torsion.rot_sym1, angle_increment) % min(360/torsion.rot_sym1, angle_increment) != 0):
                if (torsion.rot_sym2 == 1 or max(360/torsion.rot_sym2, angle_increment) % min(360/torsion.rot_sym2, angle_increment) != 0):
                    for j in range(1, 360//angle_increment):
                        torsion.rot_angles.append(j*angle_increment)
                    torsion_done[i] = 1
                    continue
            if (torsion.sym_rot_atoms1 == torsion.rot_atoms1):
                if (torsion.rot_sym1 > 1):
                    if (max(360/torsion.rot_sym1, angle_increment) % min(360/torsion.rot_sym1, angle_increment) == 0):
                        j = 1
                        while (j*angle_increment < 360/torsion.rot_sym1):
                            torsion.rot_angles.append(j*angle_increment)
                            j += 1
                        torsion_done[i] = 1
                        torsion_restricted[i] = 1
                        continue
            if (torsion.sym_rot_atoms2 == torsion.rot_atoms2):
                if (torsion.rot_sym2 > 1):
                    if (max(360/torsion.rot_sym2, angle_increment) % min(360/torsion.rot_sym2, angle_increment) == 0):
                        j = 1
                        while (j*angle_increment < 360/torsion.rot_sym2):
                            torsion.rot_angles.append(j*angle_increment)
                            j += 1
                        torsion_done[i] = 1
                        torsion_restricted[i] = 1
                        continue
            """
        
        for i, torsion1 in enumerate(torsions):
            if (torsion_done[i]):
                continue
            for j, torsion2 in enumerate(torsions[i+1:], start=i+1):
                if (torsion1.rot_sym1 > 1):
                    if (max(360/torsion1.rot_sym1, angle_increment) % min(360/torsion1.rot_sym1, angle_increment) == 0):
                        if (torsion1.sym_rot_atoms1 == torsion2.sym_rot_atoms1 or torsion1.sym_rot_atoms1 == torsion2.sym_rot_atoms2):
                            if torsion_done[j]:
                                for k in range(1, 360//angle_increment):
                                    torsion1.rot_angles.append(k*angle_increment)
                                """
                                if torsion_restricted[j]:
                                    for k in range(1, 360//angle_increment):
                                        torsion1.rot_angles.append(k*angle_increment)
                                else:
                                    k = 1
                                    while (k*angle_increment < 360/torsion1.rot_sym1):
                                        torsion1.rot_angles.append(k*angle_increment)
                                        k += 1
                                    torsion_restricted[i] = 1
                                """
                            else:
                                if len(torsion1.torsion_atoms) > len(torsion2.torsion_atoms):
                                    k = 1
                                    while (k*angle_increment < 360/torsion1.rot_sym1):
                                        torsion1.rot_angles.append(k*angle_increment)
                                        k += 1
                                    for k in range(1, 360//angle_increment):
                                        torsion2.rot_angles.append(k*angle_increment)
                                else:
                                    k = 1
                                    while (k*angle_increment < 360/torsion1.rot_sym1):
                                        torsion2.rot_angles.append(k*angle_increment)
                                        k += 1
                                    for k in range(1, 360//angle_increment):
                                        torsion1.rot_angles.append(k*angle_increment)
                            torsion_done[i] = 1
                            torsion_done[j] = 1
                            continue
                if (torsion1.rot_sym2 > 1):
                    if (max(360/torsion1.rot_sym2, angle_increment) % min(360/torsion1.rot_sym2, angle_increment) == 0):
                        if (torsion1.sym_rot_atoms2 == torsion2.sym_rot_atoms1 or torsion1.sym_rot_atoms2 == torsion2.sym_rot_atoms2):
                            if torsion_done[j]:
                                for k in range(1, 360//angle_increment):
                                        torsion1.rot_angles.append(k*angle_increment)
                                """
                                if torsion_restricted[j]:
                                    for k in range(1, 360//angle_increment):
                                        torsion1.rot_angles.append(k*angle_increment)
                                else:
                                    k = 1
                                    while (k*angle_increment < 360/torsion1.rot_sym2):
                                        torsion1.rot_angles.append(k*angle_increment)
                                        k += 1
                                    torsion_restricted[i] = 1
                                """
                            else:
                                if len(torsion1.torsion_atoms) > len(torsion2.torsion_atoms):
                                    k = 1
                                    while (k*angle_increment < 360/torsion1.rot_sym2):
                                        torsion1.rot_angles.append(k*angle_increment)
                                        k += 1
                                    for k in range(1, 360//angle_increment):
                                        torsion2.rot_angles.append(k*angle_increment)
                                else:
                                    k = 1
                                    while (k*angle_increment < 360/torsion1.rot_sym2):
                                        torsion2.rot_angles.append(k*angle_increment)
                                        k += 1
                                    for k in range(1, 360//angle_increment):
                                        torsion1.rot_angles.append(k*angle_increment)
                            torsion_done[i] = 1
                            torsion_done[j] = 1
                            continue

        """
        for i, torsion1 in enumerate(self.torsions):
            if (torsion_done[i]):
                continue
            if (torsion1.rot_sym1 == 1 or max(360/torsion1.rot_sym1, angle_increment) % min(360/torsion1.rot_sym1, angle_increment) != 0):
                if (torsion1.rot_sym2 == 1 or max(360/torsion1.rot_sym2, angle_increment) % min(360/torsion1.rot_sym2, angle_increment) != 0):
                    for j in range(1, 360//angle_increment):
                        torsion1.rot_angles.append(j*angle_increment)
                    torsion_done[i] = 1
                    continue
            if (torsion1.sym_rot_atoms1 == torsion1.rot_atoms1):
                if (torsion1.rot_sym1 > 1):
                    if (max(360/torsion.rot_sym1, angle_increment) % min(360/torsion.rot_sym1, angle_increment) == 0):
                        j = 1
                        while (j*angle_increment < 360/torsion1.rot_sym1):
                            torsion1.rot_angles.append(j*angle_increment)
                            j += 1
                        torsion_done[i] = 1
                        continue
            if (torsion1.sym_rot_atoms2 == torsion1.rot_atoms2):
                if (torsion1.rot_sym2 > 1):
                    if (max(360/torsion.rot_sym2, angle_increment) % min(360/torsion.rot_sym2, angle_increment) == 0):
                        j = 1
                        while (j*angle_increment < 360/torsion1.rot_sym2):
                            torsion1.rot_angles.append(j*angle_increment)
                            j += 1
                        torsion_done[i] = 1
                        continue
            for j, torsion2 in enumerate(self.torsions[i+1:], start=i+1):
                if (torsion_done[j]):
                    continue
                if (torsion1.rot_sym1 > 1):
                    if (max(360/torsion.rot_sym1, angle_increment) % min(360/torsion.rot_sym1, angle_increment) == 0):
                        if (torsion1.sym_rot_atoms1 == torsion2.sym_rot_atoms1 or torsion1.sym_rot_atoms1 == torsion2.sym_rot_atoms2):
                            k = 1
                            while (k*angle_increment < 360/torsion1.rot_sym1):
                                torsion2.rot_angles.append(k*angle_increment)
                                k += 1
                            torsion_done[j] = 1
                            for k in range(1, 360//angle_increment):
                                torsion1.rot_angles.append(k*angle_increment)
                            torsion_done[i] = 1
                            continue
                if (torsion1.rot_sym2 > 1):
                    if (max(360/torsion.rot_sym2, angle_increment) % min(360/torsion.rot_sym2, angle_increment) == 0):
                        if (torsion1.sym_rot_atoms2 == torsion2.sym_rot_atoms1 or torsion1.sym_rot_atoms2 == torsion2.sym_rot_atoms2):
                            k = 1
                            while (k*angle_increment < 360/torsion1.rot_sym2):
                                torsion2.rot_angles.append(k*angle_increment)
                                k += 1
                            torsion_done[j] = 1
                            for k in range(1, 360//angle_increment):
                                torsion1.rot_angles.append(k*angle_increment)
                            torsion_done[i] = 1
                            continue
        """
        """
        print()
        print(f"Inkrement: {angle_increment}")
        for torsion in torsions:
            print(f"{torsion.atom1} {torsion.atom2}")
            for angle in torsion.rot_angles:
                print(angle, end=" ")
            print()
        print()
        """
    


    def _find_torsion_group(self, connectivity: dict, torsions: list, atom: str, last_atom: str, status: dict, torsion_atoms: list) -> None:
        if status[atom] == "SEEN":
            return
        else:
            status[atom] = "SEEN"
            for neighbor in connectivity[atom]:
                if neighbor != last_atom:
                    if not atom_in_torsions(torsions, neighbor):
                        torsion_atoms.append(neighbor)
                        self._find_torsion_group(connectivity, torsions, neighbor, atom, status, torsion_atoms)
        return

    

    def rot_sym_along_bond(self, mol: Structure, rot_axis: RotationAxis, rot_atoms: list, rot_order: int) -> bool:
        for atom_i in rot_atoms:
            coords_i = mol.coords[atom_i]
            best_distance = 1.0
            symmetric_coords = rot_axis.rotate(coords_i, 360.0/float(rot_order))
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
        for rot_order in self._possible_rot_orders:
            if (self.rot_sym_along_bond(mol, rot_axis, rot_atoms, rot_order)):
                return order
        return 1