import os
import queue
import subprocess
import threading

import numpy as np
import matplotlib.pyplot as plt

from SUPRAConformer.structure import Structure
from utils.optimizer import Optimizer
from utils.analyzer import Analyzer
from utils.symmetry import Symmetry
from utils.rotationaxis import RotationAxis
from utils.bond import Bond
from utils.helper import covalence_radii_single, covalence_radii_double, get_element, \
                         increment_combinations, valences, atom_in_torsions, get_number



class ConformerGenerator:

    def __init__(self):
        self.output_folder_name = "SUPRA_Output"
        self.workdir_name = "optdir"
        self.opt_struc_name = "opt_struc.xyz"
        self.analyzer = Analyzer()
        self.optimizer = Optimizer()
        self.torsions = []
        self.central_torsions = []
        self.terminal_torsions = []
        self.methyl_torsions = []
        self.angle_increments = list()
        self.angles = list()


    def generate_conformers(
        self,
        structure: Structure,
        increment: int = 120,
        ignore_methyl: bool = False,
        ignore_terminal: bool = False,
        ignore_peptide: bool = False
    ) -> None:
        self._get_torsions(structure.bonds, structure.bond_partners)
        self._find_cycles(structure.bond_partners)
        ###
        #self._find_peptidebonds(structure.coords, structure.bond_partners)
        #self._selection_menu()
        if ignore_peptide:
            self._find_peptidebonds(structure.coords, structure.bond_partners)
        self.torsions = self.central_torsions
        if not ignore_methyl:
            self.torsions = self.torsions + self.methyl_torsions
        if not ignore_terminal:
            self.torsions = self.torsions + self.terminal_torsions
        self.angle_increments = increment_combinations[increment]
        possible_number_of_conformers = 0
        for angle in self.angle_increments:
            possible_number_of_conformers += int(np.power((360 / angle), len(self.torsions)))
        confirm = False
        while not confirm:
            user_input = input(
                f"Up to {possible_number_of_conformers} structures will be generated. Start calculation (1) or exit SUPRA (2)?: "
            )
            if user_input == "1":
                confirm = True
            elif user_input == "2":
                return
            else:
                print("Invalid input.")
        ###
        self._generation_setup(list(structure.coords.keys()), structure.bond_partners)
        print("Performing generation of conformer structures...")
        number_conformers = 0
        for increment in self.angle_increments:
            self.angles = [n * increment for n in range(int(np.round(360 / increment)))]
            self.check_rot_sym_of_torsions(structure, increment)
            number_conformers = self._combinations(
                bond_partners=structure.bond_partners, new_coords=structure.coords, counter=number_conformers
            )
        print("Generation of conformer structures done.")
        if number_conformers:
            if not os.path.exists(self.output_folder_name):
                os.makedirs(self.output_folder_name)
            self.output_folder_name = os.path.abspath(self.output_folder_name)
            for i in range(number_conformers):
                os.system(f"mv conformer{i}.xyz {self.output_folder_name}")
        return number_conformers
    


    def _get_torsions(self, bonds: list, bond_partners: dict) -> None:
        for bond in bonds:
            if bond.bond_order != 1:
                continue
            valence1 = valences[get_element(bond.atom1)]
            valence2 = valences[get_element(bond.atom2)]
            if valence1 == 1 or valence2 == 1:
                continue
            bond_partners1 = sorted(bond_partners[bond.atom1], key=lambda x: valences[get_element(x)], reverse=True)
            sum1 = sum([valences[get_element(atom)] for atom in bond_partners1[1:]])
            bond_partners2 = sorted(bond_partners[bond.atom2], key=lambda x: valences[get_element(x)], reverse=True)
            sum2 = sum([valences[get_element(atom)] for atom in bond_partners2[1:]])
            if sum1 == valence1-1:
                if get_element(bond.atom1) == "C" and \
                  [get_element(atom) for atom in bond_partners1[1:]] == ["H", "H", "H"]:
                    self.methyl_torsions.append(bond)
                else:
                    self.terminal_torsions.append(bond)
            elif sum2 == valence2-1:
                if get_element(bond.atom2) == "C" and \
                  [get_element(atom) for atom in bond_partners2[1:]] == ["H", "H", "H"]:
                    self.methyl_torsions.append(bond)
                else:
                    self.terminal_torsions.append(bond)
            else:
                self.central_torsions.append(bond)



    def _find_cycles(self, bond_partners: dict) -> None:
        atoms = list(bond_partners.keys())
        start = atoms[0]
        status = {atom: "UNKNOWN" for atom in atoms}
        self._cycle_detection(bond_partners, start, start, status)



    # find rings in structure
    # based on Depth-First-Search (DFS) with coloring method
    # UNKNOWN: node not visited yet
    # VISITED: node visited but not done yet
    # KNOWN: node done
    def _cycle_detection(self, bond_partners: dict, start: str, last: str, status: dict, ancestors: dict = {}) -> None:
        # base case 1: node done or terminal atom
        if status[start] == "KNOWN" or get_element(start) in ["H", "F", "Cl", "Br", "I"]:
            return
        # base case 2: ring found (visited node again on a path where it has already been visited)
        elif status[start] == "VISITED":
            # loop over all bond in list of rotatable bonds
            for torsion in self.central_torsions:
                # if bond is identical with first bond of current cycle delete it from list of
                # rotatable bonds
                if (torsion.atom1 == start and torsion.atom2 == last) or (torsion.atom2 == start and torsion.atom1 == last):
                    self.central_torsions.remove(torsion)
                    # if bond has been found in list break loop
                    break
            cur = last
            # get all atoms of current ring by backtracking
            # loop while full ring has not been traversed yet
            while cur != start:
                # get current atom in ring
                before = cur
                # get next atom in ring
                cur = ancestors[cur]
                # loop over all rotatable bonds
                for torsion in self.central_torsions:
                    # if current and next form bond which is in list of rotatable bond delete it from the list
                    if (torsion.atom1 == before and torsion.atom2 == cur) or (torsion.atom2 == before and torsion.atom1 == cur):
                        self.central_torsions.remove(torsion)
                        # if bond has been found in list break loop
                        break
        # case: node not visited yet
        else:
            # store last visited node
            ancestors[start] = last
            # mark current node as visited
            status[start] = "VISITED"
            # visit all neighbouring nodes of current node in sence of a dfs
            for bond_partner in bond_partners[start]:
                # only traverse through chain, ignore terminal atoms
                if not get_element(bond_partner) in ["H", "F", "Cl", "Br", "I"]:
                    # no jumping forth and back between current and last visited node
                    if not bond_partner == ancestors[start]:
                        self._cycle_detection(bond_partners, bond_partner, start, status, ancestors)
            status[start] = "KNOWN"



    # delete peptide bonds from list of rotatable bonds (torsions)
    def _find_peptidebonds(self, coords: dict, bond_partners: dict) -> None:
        peptidebonds = []
        double_bond = covalence_radii_double["C"] + covalence_radii_double["O"] + 0.08
        # loop over every rotatable bond
        for bond in self.central_torsions:
            # check for characteristic C=O + C-NRR' bonding environment
            # first check: C and N included in bond?
            if sorted([get_element(bond.atom1), get_element(bond.atom2)]) == ["C", "N"]:
                # second check: C with bonding partner O?
                if get_element(bond.atom1) == "C":
                    C = bond.atom1
                else:
                    C = bond.atom2
                for atom3 in bond_partners[C]:
                    element3 = get_element(atom3)
                    if element3 == "O":
                        # third check: double bond between C and O?
                        p1 = coords[C]
                        p2 = coords[atom3]
                        distance = np.linalg.norm(p1 - p2)
                        if distance <= double_bond:
                            # found peptide bond, remove bond from list of rotatable bonds
                            peptidebonds.append(bond)
                            break
        for bond in peptidebonds:
            self.central_torsions.remove(bond)
        #if peptidebonds:
        #    while True:
        #        rotatable = input("Consider peptide bonds as rotatable (1) or not (2)?: ")
        #        if rotatable == "1":
        #            break
        #        elif rotatable == "2":
        #            for bond in peptidebonds:
        #                self.central_torsions.remove(bond)
        #            break
        #        else:
        #            print("Invalid input.")
    


    def _selection_menu(self):
        start_calculation = False
        while not start_calculation:
            while True:
                mode_input = input(
                    "Consider rotatable all bonds to terminal groups like "
                    "-CH3, -NH2, -OH (1) "
                    "or ignore them (2) "
                    "or ignore just ignore bonds to methyl groups (3)? "
                )
                if mode_input == "1":
                    self.torsions = self.central_torsions + self.terminal_torsions + self.methyl_torsions
                    break
                elif mode_input == "2":
                    self.torsions = self.central_torsions
                    break
                elif mode_input == "3":
                    self.torsions = self.central_torsions + self.terminal_torsions
                    break
                else:
                    print("Invalid input.")
            confirm_increment = False
            while not confirm_increment:
                increment_input = input(f"Type in an angle increment (30, 45, 60, 90, 120 or 180 in degrees): ")
                if increment_input in ["30", "45", "60", "90", "120", "180"]:
                    self.angle_increments = increment_combinations[int(increment_input)]
                    possible_number_of_conformers = 0
                    for angle in self.angle_increments:
                        possible_number_of_conformers += int(np.power((360 / angle), len(self.torsions)))
                    confirm = -1
                    while True:
                        confirm = input(f"Up to {possible_number_of_conformers} will be generated. Start calculation (1) or restart selection (2)?: ")
                        if confirm == "1":
                            confirm_increment = True
                            start_calculation = True
                            break
                        elif confirm == "2":
                            confirm_increment = True
                            break
                        else:
                            print("Invalid input.")
                else:
                    print("Invalid input.")



    # graphical output of the input molecule structure with numbering of the atoms and coloring of bonds
    def show_structure(self, coords: dict, bonds: list, selection_thread: threading.Thread = None) -> None:
        # lists for coordinates, bonds and labels of chain atoms
        central_xcoords = list()
        central_ycoords = list()
        central_zcoords = list()
        central_labels = list()
        # lists for coordinates, bonds and labels of terminal atoms
        terminal_xcoords = list()
        terminal_ycoords = list()
        terminal_zcoords = list()
        terminal_labels = list()
        # lists for plots
        terminal_torsions = list()
        central_torsions = list()
        no_torsions_central = list()
        no_torsions_terminal = list()
        # loop over all atoms
        for atom1, coords1 in coords.items():
            # get coordinates of first atom of a bond
            x1, y1, z1 = coords1
            # if atom is terminal raise flag
            terminal1 = get_element(atom1) in ["H", "F", "Cl", "Br", "I"]
            # if flag is raised store coordinates and label of atom in according lists of terminal atoms
            if terminal1:
                # store coordinates
                terminal_xcoords.append(x1)
                terminal_ycoords.append(y1)
                terminal_zcoords.append(z1)
                # store label
                terminal_labels.append(atom1)
            # if flag is not raised save coordinates and label of atom in according lists of chain atoms
            else:
                # store coordinates
                central_xcoords.append(x1)
                central_ycoords.append(y1)
                central_zcoords.append(z1)
                # store label
                central_labels.append(atom1)
        for bond in bonds:
            atom1, atom2 = bond
            x1, y1, z1 = coords[atom1]
            x2, y2, z2 = coords[atom2]
            bond_coords = [[x1, x2], [y1, y2], [z1, z2]]
            if bond in self.terminal_torsions:
                terminal_torsions.append(bond_coords)
            elif bond in self.central_torsions:
                central_torsions.append(bond_coords)
            elif get_element(atom1) in ["H", "F", "Cl", "Br", "I"] or \
                 get_element(atom2) in ["H", "F", "Cl", "Br", "I"]:
                no_torsions_terminal.append(bond_coords)
            else:
                no_torsions_central.append(bond_coords)
        # create new empty 3D plot
        ax = plt.axes(projection="3d")
        # hide axis, just molecule structure should be shown
        ax.set_axis_off()
        # plot chain atoms at coordinates of chain atoms (dots, black)
        ax.scatter(central_xcoords, central_ycoords, central_zcoords, color="black")
        # plot terminal atoms at coordinates of terminal atoms (dots, light gray)
        ax.scatter(terminal_xcoords, terminal_ycoords, terminal_zcoords, color="lightgrey")
        # plot labels of chain atoms next to atoms (black)
        for iterator, atom in enumerate(central_labels):
            ax.text(central_xcoords[iterator], central_ycoords[iterator], central_zcoords[iterator], atom)
        # plot labels of terminal atoms next to atoms (light gray)
        for iterator, atom in enumerate(terminal_labels):
            ax.text(terminal_xcoords[iterator], terminal_ycoords[iterator], terminal_zcoords[iterator], atom,
                    color="lightgrey")
        # declaration of indices for better readability of code
        atom1 = 0
        atom2 = 1
        x = 0
        y = 1
        z = 2
        # plot rotatable bonds to terminal groups in blue
        for bond in terminal_torsions:
            ax.plot([bond[x][atom1], bond[x][atom2]], [bond[y][atom1], bond[y][atom2]],
                    [bond[z][atom1], bond[z][atom2]], color="blue")
        # plot rotatable, central bonds in red
        for bond in central_torsions:
            ax.plot([bond[x][atom1], bond[x][atom2]], [bond[y][atom1], bond[y][atom2]],
                    [bond[z][atom1], bond[z][atom2]], color="red")
        # plot non-rotatable bonds in black
        for bond in no_torsions_central:
            ax.plot([bond[x][atom1], bond[x][atom2]], [bond[y][atom1], bond[y][atom2]],
                    [bond[z][atom1], bond[z][atom2]], color="black")
        # plot terminal bonds in light gray
        for bond in no_torsions_terminal:
            ax.plot([bond[x][atom1], bond[x][atom2]], [bond[y][atom1], bond[y][atom2]],
                    [bond[z][atom1], bond[z][atom2]], color="lightgrey")
        # if there is a parallel thread with bond selection special handling is necessary
        if selection_thread:
            plt.draw()
            # show plot while parallel thread with bond selection is alive
            while selection_thread.is_alive():
                plt.pause(1)
        else:
            plt.show()



    # zählen der jeweils zu rotierenden Atome auf beiden Seiten der Bindung, Hinzufügen zu torsion_atoms von der
    # Seite welche weniger Atome enthält
    def _generation_setup(self, atoms: list, bond_partners: dict) -> None:
        # Durchführung für jede gefundene und gefilterte bzw. ausgewählte Rotationsbindung
        for torsion in self.torsions:
            # zählen der linken Seite (entspricht erstgenanntem Atom der Bindung in torsions)
            status = {atom: "UNKNOWN" for atom in atoms}
            self._torsion_atom_counter(bond_partners, torsion.atom1, torsion.atom2, status, torsion.rot_atoms1)
            # zählen der linken Seite (entspricht letztgenanntem Atom der Bindung in torsions)
            status = {atom: "UNKNOWN" for atom in atoms}
            self._torsion_atom_counter(bond_partners, torsion.atom2, torsion.atom1, status, torsion.rot_atoms2)
            # wenn links weniger oder genauso viele Atome wie rechts wird links gedreht, d.h. Atome der linken
            # Seite zu torsions_atoms hinzufügen
            if len(torsion.rot_atoms1) <= len(torsion.rot_atoms2):
                torsion.torsion_atoms = torsion.rot_atoms1
            # wenn rechts weniger  Atome wie links wird links gedreht, d.h. Atome der rechten Seite zu
            # torsions_atoms hinzufügen
            else:
                torsion.torsion_atoms = torsion.rot_atoms2



    # Hilfsfunktion, wird von rotation_setup aufgerufen und zählt zu rotierende Atome ausgehend vom Atom einer
    # Bindung bei Rotation um diese
    #
    # basiert auf Depth-First-Search (DFS) auf ungerichtetem Graphen mit Färbemethode (analog zu find_cycles)
    # UNKNOWN: Knoten noch nicht besucht
    # SEEN: Knoten besucht
    def _torsion_atom_counter(
            self, bond_partners: dict, atom: str, last_atom: str, status: dict, torsion_atoms: list
    ) -> list:
        # Atom bereits betrachtet
        if status[atom] == "SEEN":
            return
        # Atom noch nicht betrachtet und nicht Ende einer Kette
        else:
            # Atom als betrachtet markieren
            status[atom] = "SEEN"
            # rekursiver Aufruf der Funktion für alle Bindungspartner von Atom, solange diesem nicht dem
            # Atom des vorherigen Funktionsaurufs entsprechen
            for neighbor in bond_partners[atom]:
                if neighbor != last_atom:
                    torsion_atoms.append(neighbor)
                    self._torsion_atom_counter(
                        bond_partners, neighbor, atom, status, torsion_atoms
                    )
        # Anzahl der zu rotierenden Atome der entsprechenden Seite der Bindung zurückgeben
        return
    


    def check_rot_sym_of_torsions(self, mol: Structure, angle_increment: int) -> None:
        sym = Symmetry()

        torsion_done = [0 for _ in self.torsions]
        torsion_restricted = [0 for _ in self.torsions]

        for i, torsion in enumerate(self.torsions):
            atom1 = torsion.atom1
            atom2 = torsion.atom2
            torsion.rot_angles = [0]
            # rotational symmetry left side of torsion bond
            status = {atom: "UNKNOWN" for atom in mol.coords.keys()}
            self._get_torsion_group(mol.bond_partners, atom1, atom2, status, torsion.sym_rot_atoms1)
            torsion.sym_rot_atoms1 = sorted(torsion.sym_rot_atoms1, key=lambda label: get_number(label))
            torsion.rot_sym1 = sym.rot_order_along_bond(mol, torsion.sym_rot_atoms1, mol.coords[atom1], mol.coords[atom2])
            # rotational symmetry right side of torsion bond
            status = {atom: "UNKNOWN" for atom in mol.coords.keys()}
            self._get_torsion_group(mol.bond_partners, torsion.atom2, torsion.atom1, status, torsion.sym_rot_atoms2)
            torsion.sym_rot_atoms2 = sorted(torsion.sym_rot_atoms2, key=lambda label: get_number(label))
            torsion.rot_sym2 = sym.rot_order_along_bond(mol, torsion.sym_rot_atoms2, mol.coords[atom1], mol.coords[atom2])
            # assign rotation angles if already possible
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
                            torsion1.rot_angles.append(j*angle_increment)
                            j += 1
                        torsion_done[i] = 1
                        torsion_restricted[i] = 1
                        continue
            if (torsion.sym_rot_atoms2 == torsion.rot_atoms2):
                if (torsion.rot_sym2 > 1):
                    if (max(360/torsion.rot_sym2, angle_increment) % min(360/torsion.rot_sym2, angle_increment) == 0):
                        j = 1
                        while (j*angle_increment < 360/torsion.rot_sym2):
                            torsion1.rot_angles.append(j*angle_increment)
                            j += 1
                        torsion_done[i] = 1
                        torsion_restricted[i] = 1
                        continue
        
        for i, torsion1 in enumerate(self.torsions):
            if (torsion_done[i]):
                continue
            for j, torsion2 in enumerate(self.torsions[i+1:], start=i+1):
                if (torsion1.rot_sym1 > 1):
                    if (max(360/torsion.rot_sym1, angle_increment) % min(360/torsion.rot_sym1, angle_increment) == 0):
                        if (torsion1.sym_rot_atoms1 == torsion2.sym_rot_atoms1 or torsion1.sym_rot_atoms1 == torsion2.sym_rot_atoms2):
                            if torsion_done[j]:
                                if torsion_restricted[j]:
                                    for k in range(1, 360//angle_increment):
                                        torsion1.rot_angles.append(k*angle_increment)
                                else:
                                    k = 1
                                    while (k*angle_increment < 360/torsion1.rot_sym1):
                                        torsion1.rot_angles.append(k*angle_increment)
                                        k += 1
                                    torsion_restricted[i] = 1
                            else:
                                if len(torsion1.torsion_atoms) > len(torsion.torsion_atoms):
                                    k = 1
                                    while (k*angle_increment < 360/torsion1.rot_sym1):
                                        torsion1.rot_angles.append(k*angle_increment)
                                        k += 1
                                    for k in range(1, 360//angle_increment):
                                        torsion2.rot_angles.append(k*angle_increment)
                                    torsion_restricted[i] = 1
                                else:
                                    k = 1
                                    while (k*angle_increment < 360/torsion1.rot_sym1):
                                        torsion2.rot_angles.append(k*angle_increment)
                                        k += 1
                                    for k in range(1, 360//angle_increment):
                                        torsion1.rot_angles.append(k*angle_increment)
                                    torsion_restricted[j] = 1
                            torsion_done[i] = 1
                            torsion_done[j] = 1
                            continue
                if (torsion1.rot_sym2 > 1):
                    if (max(360/torsion.rot_sym2, angle_increment) % min(360/torsion.rot_sym2, angle_increment) == 0):
                        if (torsion1.sym_rot_atoms2 == torsion2.sym_rot_atoms1 or torsion1.sym_rot_atoms2 == torsion2.sym_rot_atoms2):
                            if torsion_done[j]:
                                if torsion_restricted[j]:
                                    for k in range(1, 360//angle_increment):
                                        torsion1.rot_angles.append(k*angle_increment)
                                else:
                                    k = 1
                                    while (k*angle_increment < 360/torsion1.rot_sym2):
                                        torsion1.rot_angles.append(k*angle_increment)
                                        k += 1
                                    torsion_restricted[i] = 1
                            else:
                                if len(torsion1.torsion_atoms) > len(torsion.torsion_atoms):
                                    k = 1
                                    while (k*angle_increment < 360/torsion1.rot_sym2):
                                        torsion1.rot_angles.append(k*angle_increment)
                                        k += 1
                                    for k in range(1, 360//angle_increment):
                                        torsion2.rot_angles.append(k*angle_increment)
                                    torsion_restricted[i] = 1
                                else:
                                    k = 1
                                    while (k*angle_increment < 360/torsion1.rot_sym2):
                                        torsion2.rot_angles.append(k*angle_increment)
                                        k += 1
                                    for k in range(1, 360//angle_increment):
                                        torsion1.rot_angles.append(k*angle_increment)
                                    torsion_restricted[j] = 1
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
        for i, torsion in enumerate(self.torsions):
            print(f"{torsion.atom1} {torsion.atom2}")
            for angle in torsion.rot_angles:
                print(angle, end=" ")
            print()
        print()
        """



    def _get_torsion_group(self, connectivity: dict, atom: str, last_atom: str, status: dict, torsion_atoms: list):
        if status[atom] == "SEEN":
            return
        else:
            status[atom] = "SEEN"
            for neighbor in connectivity[atom]:
                if neighbor != last_atom:
                    if not atom_in_torsions(self.torsions, neighbor):
                        torsion_atoms.append(neighbor)
                        self._get_torsion_group(connectivity, neighbor, atom, status, torsion_atoms)
        return



    # calculates all possible conformer structures and generates an output file for every conformer structure without
    # internal clash
    # the output file is an input file for a geometry optimization with ORCA
    def _combinations(self, bond_partners: dict, new_coords: dict, counter: int, index: int = 0) -> int:
        # base case, new torsion angle for every angle has been calculated
        if index == len(self.torsions):
            # sofern keine strukturinternen Clashes hinzufügen zur Liste erfolgreich erzeugter Konformerstrukturen
            # check new structure for internal clashes
            if not self._clashes(bond_partners, new_coords):
                self.optimizer.uff_structure_optimization(new_coords, counter)
                return counter+1
            else:
                return counter
        # es wurden noch nicht alle Torsionswinkel berechnet, Bindung index in torsions ist an der Reihe
        else:
            # Punkte initialisieren, welche die Rotationsachse definieren
            axis_vec1 = new_coords[self.torsions[index].atom1]
            axis_vec2 = new_coords[self.torsions[index].atom2]
            # jeden möglichen Torsionswinkel für Bindung durchgehen
            #for angle in self.angles:
            for angle in self.torsions[index].rot_angles:
                new_coords_copy = new_coords.copy()
                for atom in self.torsions[index].torsion_atoms:
                    new_coords_copy[atom] = RotationAxis.rotate_atom(axis_vec1, axis_vec2, new_coords[atom], angle)
                # rekursiver Aufruf für nächste Bindung
                counter = self._combinations(bond_partners, new_coords_copy, counter, index+1)
            return counter



    def _clashes(self, bond_partners: dict, new_coords: dict) -> bool:
        for atom1, coords1 in new_coords.items():
            element1 = get_element(atom1)
            for atom2, coords2 in new_coords.items():
                if not atom1 == atom2:
                    element2 = get_element(atom2)
                    distance = np.linalg.norm(coords1 - coords2)
                    min_distance = covalence_radii_single[element1] + covalence_radii_single[element2] + 0.15#0.08
                    if distance < min_distance:
                        if self._distant(bond_partners, atom1, atom2):
                            return True
        return False



    def _distant(self, bond_partners: dict, atom1: str, atom2: str) -> bool:
        neighbors = queue.Queue()
        distance = 1
        for neighbor in bond_partners[atom1]:
            neighbors.put(neighbor)
        while distance < 4 and neighbors:
            num_neighbors = neighbors.qsize()
            while num_neighbors:
                curr = neighbors.get()
                if curr == atom2:
                    return False
                else:
                    for neighbor in bond_partners[curr]:
                        neighbors.put(neighbor)
                num_neighbors -= 1
            distance += 1
        return True



    def output_coords(self, coords: dict, counter: int) -> None:
        output_folder = "../Output/"
        if not os.path.exists(output_folder):
            os.makedirs(output_folder)
        with open(f"{output_folder}conformer{counter}.xyz", "w") as outfile:
            for atom, (x, y, z) in coords.items():
                element = get_element(atom)
                print(f"{element}\t{x:18.15f}\t{y:18.15f}\t{z:18.15f}", file=outfile)



    def output_opt(self, coords: dict, counter: int) -> None:
        output_folder = "../Output/"
        if not os.path.exists(output_folder):
            os.makedirs(output_folder)
        with open(f"{output_folder}conformer{counter}", "w") as outfile:
            outfile.write("! BP86 SVP D4 TightSCF TightOpt NumFreq PAL4\n\n")
            outfile.write("*xyz 0 1\n")
            print("! BP86 SVP D4 TightSCF TightOpt NumFreq PAL4", file=outfile, end="\n\n")
            print("*xyz 0 1", file=outfile)
            for atom, (x, y, z) in coords:
                print(f"{get_element(atom)}\t{x:18.15f}\t{y:18.15f}\t{z:18.15f}", file=outfile)
            print("*", file=outfile)
