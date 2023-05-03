import os
import queue
import subprocess
import threading

import numpy as np
import matplotlib.pyplot as plt

from SUPRAConformer.structure import Structure
from utils.analyzer import Analyzer
from utils.helper import covalence_radii_single, covalence_radii_double, get_element, increment_combinations, valences



class ConformerGenerator:

    def __init__(self):
        self.output_folder_name = "SUPRA_Output"
        self.workdir_name = "optdir"
        self.opt_struc_name = "opt_struc.xyz"
        self.analyzer = Analyzer()
        self.torsions = []
        self.central_torsions = []
        self.terminal_torsions = []
        self.methyl_torsions = []
        self.coords = {}
        # atoms to be rotated in rotation around bond at same index in list of rotatable bonds
        # torsion_atoms = [[E0, E1, E2, ... ], [E3, E4, ... ], ... ]
        self.torsion_atoms = []
        # increment by which torsion angles are varied for conformer generation (default = 0)
        self.angle_increments = list()
        self.angles = list()


    def generate_conformers(self, structure: Structure) -> None:
        self._get_torsions(structure.bonds, structure.bond_partners, structure.bond_orders)
        self._find_cycles(structure.bond_partners)
        self._find_peptidebonds(structure.coords, structure.bond_partners)
        ### UNDER CONSTRUCTION ###
        #selection = self._selection_menu(structure)
        # stattdessen derzeit das:
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
        start_calculation = False
        while not start_calculation:
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
                        start_calculation = True
                        break
                    elif confirm == "2":
                        break
                    else:
                        print("Invalid input.")
            else:
                print("Invalid input.")
        ##########################
        self._generation_setup(list(structure.coords.keys()), structure.bond_partners)
        number_conformers = 0
        for increment in self.angle_increments:
            self.angles = [n * increment for n in range(int(np.round(360 / increment)))]
            number_conformers = self._combinations(
                bond_partners=structure.bond_partners, new_coords=structure.coords, counter=number_conformers
            )
        if number_conformers:
            if not os.path.exists(self.output_folder_name):
                os.makedirs(self.output_folder_name)
            self.output_folder_name = os.path.abspath(self.output_folder_name)
            for i in range(number_conformers):
                folder = os.path.abspath(f"{self.workdir_name}{i}")
                opt_struc = os.path.abspath(os.path.join(folder, self.opt_struc_name))
                #new_conformer = Structure(opt_struc)
                #double_flag = False
                #for conformer_file in os.listdir(self.output_folder_name):
                #    conformer_file = os.path.abspath(os.path.join(self.output_folder_name, conformer_file))    
                #    conformer = Structure(conformer_file)
                #    if self.analyzer.doubles(new_conformer, conformer):
                #        double_flag = True
                #        break
                #if not double_flag:
                new_conformer_file = os.path.join(self.output_folder_name, f"conformer{i}.xyz")
                os.system(f"mv {opt_struc} {new_conformer_file}")
                os.system(f"rm -rf {folder}")
            print(f"{len(os.listdir(self.output_folder_name))} conformers have been generated.")
        else:
            print(f"No conformers could be generated.")


    def _get_torsions(self, bonds: list, bond_partners: dict, bond_orders: dict) -> None:
        for bond in bonds:
            atom1, atom2 = bond
            valence1 = valences[get_element(atom1)]
            valence2 = valences[get_element(atom2)]
            if bond_orders[bond] != 1 or valence1 == 1 or valence2 == 1:
                continue
            else:
                bond_partners1 = sorted(bond_partners[atom1], key=lambda x: valences[get_element(x)], reverse=True)
                sum1 = sum([valences[get_element(atom)] for atom in bond_partners1[1:]])
                bond_partners2 = sorted(bond_partners[atom2], key=lambda x: valences[get_element(x)], reverse=True)
                sum2 = sum([valences[get_element(atom)] for atom in bond_partners2[1:]])
                if sum1 == valence1-1 :
                    if get_element(atom1) == "C" and \
                      [get_element(atom) for atom in bond_partners1[1:]] == ["H", "H", "H"]:
                        self.methyl_torsions.append(bond)
                    else:
                        self.terminal_torsions.append(bond)
                elif sum2 == valence2-1:
                    if get_element(atom2) == "C" and \
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
                if torsion == (start, last) or torsion == (last, start):
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
                    if torsion == (before, cur) or torsion == (cur, before):
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
        # loop over every rotatable bond
        for bond in self.central_torsions:
            atom1 = bond[0]
            atom2 = bond[1]
            element1 = get_element(atom1)
            # check for characteristic C=O + C-NHR bonding environment
            # first check: C and N included in bond?
            if sorted([atom1, atom2]) == ["C", "N"]:
                # second check: C with bonding partner O?
                for atom3 in bond_partners[atom1]:
                    element3 = get_element(atom3)
                    if element3 == "O":
                        # third check: double bond between C and O?
                        p1 = coords[atom1]
                        p2 = coords[atom3]
                        distance = np.linalg.norm(p1 - p2)
                        double_bond = covalence_radii_double[element1] + covalence_radii_double[element3] + 0.08
                        if distance <= double_bond:
                            # found peptide bond, remove bond from list of rotatable bonds
                            peptidebonds.append(bond)
                            break
        if peptidebonds:
            while True:
                rotatable = input("Consider peptide bonds as rotatable (1) or not (2)?: ")
                if rotatable == "1":
                    break
                elif rotatable == "2":
                    for bond in peptidebonds:
                        self.central_torsions.remove(bond)
                    break
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
        for bond in self.torsions:
            torsion_atoms_left = []
            torsion_atoms_right = []
            # zählen der linken Seite (entspricht erstgenanntem Atom der Bindung in torsions)
            status = {atom: "UNKNOWN" for atom in atoms}
            left = self._torsion_atom_counter(
                bond_partners, bond[0], bond[1], status, -1, torsion_atoms_left
            )
            # zählen der linken Seite (entspricht letztgenanntem Atom der Bindung in torsions)
            status = {atom: "UNKNOWN" for atom in atoms}
            right = self._torsion_atom_counter(
                bond_partners, bond[1], bond[0], status, -1, torsion_atoms_right
            )
            # wenn links weniger oder genauso viele Atome wie rechts wird links gedreht, d.h. Atome der linken
            # Seite zu torsions_atoms hinzufügen
            if left <= right:
                self.torsion_atoms.append(torsion_atoms_left)
            # wenn rechts weniger  Atome wie links wird links gedreht, d.h. Atome der rechten Seite zu
            # torsions_atoms hinzufügen
            else:
                self.torsion_atoms.append(torsion_atoms_right)


    # Hilfsfunktion, wird von rotation_setup aufgerufen und zählt zu rotierende Atome ausgehend vom Atom einer
    # Bindung bei Rotation um diese
    #
    # basiert auf Depth-First-Search (DFS) auf ungerichtetem Graphen mit Färbemethode (analog zu find_cycles)
    # UNKNOWN: Knoten noch nicht besucht
    # SEEN: Knoten besucht
    def _torsion_atom_counter(
            self, bond_partners: dict, atom: str, last_atom: str,
            status: dict, counter: int, torsion_atoms: list
    ) -> int:
        # Atom bereits betrachtet
        if status[atom] == "SEEN":
            return counter
        # Atom ist Terminalatom, d.h. Ende einer Kette
        elif get_element(atom) in ["H", "F", "Cl", "Br", "I"]:
            # Atom als betrachtet markieren
            status[atom] = "SEEN"
            # Atom zu torsion_atoms hinzufügen
            torsion_atoms.append(atom)
            return counter + 1
        # Atom noch nicht betrachtet und nicht Ende einer Kette
        else:
            # Atom als betrachtet markieren
            status[atom] = "SEEN"
            # Atom zu torsion_atoms hinzufügen
            torsion_atoms.append(atom)
            # Zähler aktualisieren
            counter += 1
            # rekursiver Aufruf der Funktion für alle Bindungspartner von Atom, solange diesem nicht dem
            # Atom des vorherigen Funktionsaurufs entsprechen
            for neighbor in bond_partners[atom]:
                if neighbor != last_atom:
                    counter = self._torsion_atom_counter(
                        bond_partners, neighbor, atom, status, counter, torsion_atoms
                    )
        # Anzahl der zu rotierenden Atome der entsprechenden Seite der Bindung zurückgeben
        return counter


    # calculates all possible conformer structures and generates an output file for every conformer structure without
    # internal clash
    # the output file is an input file for a geometry optimization with ORCA
    def _combinations(self, bond_partners: dict, new_coords: dict, counter: int, index: int = 0, optmode: str = "uff") -> int:
        # base case, new torsion angle for every angle has been calculated
        if index == len(self.torsions):
            # sofern keine strukturinternen Clashes hinzufügen zur Liste erfolgreich erzeugter Konformerstrukturen
            # check new structure for internal clashes
            if not self._clashes(bond_partners, new_coords): # wird so nicht mehr funktionieren
                workdir = f"{self.workdir_name}{counter}"
                os.makedirs(workdir)
                workdir = os.path.abspath(workdir)
                opt_struc = os.path.join(workdir, self.opt_struc_name)
                new_xyz_file = os.path.join(workdir, "struc.xyz")
                os.system(f"touch {new_xyz_file}")
                with open(new_xyz_file, "w") as optfile:
                    number_of_atoms = len(new_coords.keys())
                    print(number_of_atoms, file=optfile, end="\n\n")
                    for atom, (x, y, z) in new_coords.items():
                        print(f"{get_element(atom)}\t{x}\t{y}\t{z}", file=optfile)
                if optmode == "xtb":
                    opt_struc = os.path.join(workdir, self.opt_struc_name)
                    subprocess.run(args=["xtb", "--opt", "sloppy", new_xyz_file], cwd=workdir, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)
                    subprocess.run(args=["mv", "xtbopt.xyz", opt_struc], cwd=workdir)
                elif optmode == "uff":
                    # TODO: MAKE NUMBER OF OPTIMIZATION CYCLES DEPENDENT FROM SIZE OF MOLECULE (NUMBER OF ATOMS)
                    control_file = os.path.join(workdir, "control")
                    coord_file = os.path.join(workdir, "coord")
                    os.system(f"touch {control_file}")
                    with open(control_file, "w") as control:
                        print("$symmetry c1", file=control)
                        print("$uff", file=control)
                        print("      2500         1          0 ! maxcycle,modus,nqeq", file=control)
                        print("    111111                      ! iterm", file=control)
                        print("  0.10D-07  0.10D-04            ! econv,gconv", file=control)
                        print("      0.00  1.10                ! qtot,dfac", file=control)
                        print("  0.10D+03  0.10D-04       0.30 ! epssteep,epssearch,dqmax", file=control)
                        print("        25      0.10       0.00 ! mxls,dhls,ahls", file=control)
                        print("      1.00      0.00       0.00 ! alpha,beta,gamma", file=control)
                        print("         F         F          F ! transform,lnumhess,lmd", file=control)
                        print("$end", file=control)
                    with open(coord_file, "w") as f:
                        subprocess.run(args=["x2t", new_xyz_file, ">", coord_file], cwd=workdir, stdout=f, stderr=subprocess.DEVNULL)
                    subprocess.run(args=["uff"], cwd=workdir, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)
                    if os.path.isfile(os.path.join(workdir, "not.uffconverged")):
                        subprocess.run(args=["xtb", "--opt", "--gfnff", new_xyz_file], cwd=workdir, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)
                        os.system(f"mv {os.path.join(workdir, 'xtbopt.xyz')} {opt_struc}")
                    else:
                        #TODO: HOW TO INCLUDE ENERGY IN SECOND LINE HERE
                        with open(opt_struc, "w") as f:
                            subprocess.run(args=["t2x", coord_file, ">", opt_struc], cwd=workdir, stdout=f, stderr=subprocess.DEVNULL)
                return counter+1
            else:
                return counter
        # es wurden noch nicht alle Torsionswinkel berechnet, Bindung index in torsions ist an der Reihe
        else:
            # Punkte initialisieren, welche die Rotationsachse definieren
            axis_vec1 = new_coords[self.torsions[index][0]]
            axis_vec2 = new_coords[self.torsions[index][1]]
            # jeden möglichen Torsionswinkel für Bindung durchgehen
            for angle in self.angles:
                angle = np.deg2rad(angle)
                new_coords_copy = new_coords.copy()
                for atom in self.torsion_atoms[index]:
                    coords = new_coords[atom]
                    axis = axis_vec2 - axis_vec1
                    axis = axis / np.linalg.norm(axis)
                    coords = coords - axis_vec1
                    coords = np.dot(axis, np.dot(axis, coords)) \
                             + np.cos(angle) * np.cross(np.cross(axis, coords), axis) \
                             + np.sin(angle) * np.cross(axis, coords)
                    coords = coords + axis_vec1
                    new_coords_copy[atom] = coords
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
                    min_distance = covalence_radii_single[element1] + covalence_radii_single[element2] + 0.08
                    if distance < min_distance:
                        if self.__distant(bond_partners, atom1, atom2):
                            return True
        return False


    def __distant(self, bond_partners: dict, atom1: str, atom2: str) -> bool:
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
