from Helper import covalence_radii, get_element, get_number, rotation  # covalence radii,
                                                                       # atom label to element symbol converter,
                                                                       # atom label to number converter,
                                                                       # rotation around arbitrary axis in 3D
import matplotlib.pyplot as plt  # plot molecule in 3D-diagram
import threading  # selection menu and 3D-plot of molecule in parallel threads
import numpy as np  # sqrt, power, linalg, dot, deg2grad, cross, cos, sin, array
import os  # path, makedirs
import queue  # Queue
from Structure import Structure



class ConformerGenerator:

    def __init__(self):
        # number of rings in input structure (default = 0)
        self.cycle_number = 0

        # atoms in rings, index of list of atoms in cycles denotes number of ring
        # cycles = [[E0, E1, E3, ... ], [E0, E4, E5, ... ], ... ]
        self.cycles = list()

        # atoms to be rotated in rotation around bond at same index in list of rotatable bonds
        # torsion_atoms = [[E0, E1, E2, ... ], [E3, E4, ... ], ... ]
        self.torsion_atoms = list()

        # increment by which torsion angles are varied for conformer generation (default = 0)
        self.angle_increment = 120

        self.angles = list()

        self.conformers = list()


    def find_cycles(self, molecule: Structure):
        start = next(iter(molecule.structure))
        status = {atom: "UNKNOWN" for atom in molecule.structure}
        self._cycle_detection(molecule, start, start, status)


    # find rings in structure
    # based on Depth-First-Search (DFS) with coloring method
    # UNKNOWN: node not visited yet
    # VISITED: node visited but not done yet
    # KNOWN: node done
    def _cycle_detection(self, molecule: Structure, start: str, last: str, status: dict, ancestors: dict = dict()):
        # base case 1: node done or terminal atom
        if status[start] == "KNOWN" or get_element(start) in ["H", "F", "Cl", "Br", "I"]:
            return

        # base case 2: ring found (visited node again on a path where it has already been visited)
        elif status[start] == "VISITED":
            # create new list/storage for ring atoms of current ring in cycles
            self.cycles.append([])

            # store first atom of ring in list in cycles
            self.cycles[self.cycle_number].append(last)
            # while loop (see below starts with bond of last atom and atom before that, therefore this loop
            # to check bond of current atom and last atom)
            # loop over all bond in list of rotatable bonds
            for torsion in molecule.torsions:
                # if bond is identical with first bond of current cycle delete it from list of
                # rotatable bonds
                if torsion == [start, last] or torsion == [last, start]:
                    molecule.torsions.remove(torsion)
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
                # store current atom of ring in list in cycles
                self.cycles[self.cycle_number].append(cur)
                # loop over all rotatable bonds
                for torsion in molecule.torsions:
                    # if current and next form bond which is in list of rotatable bond delete it from the list
                    if torsion == [before, cur] or torsion == [cur, before]:
                        molecule.torsions.remove(torsion)
                        # if bond has been found in list break loop
                        break
            # increase number of rings in structure by one
            self.cycle_number += 1

        # case: node not visited yet
        else:
            # store last visited node
            ancestors[start] = last
            # mark current node as visited
            status[start] = "VISITED"

            # visit all neighbouring nodes of current node in sence of a dfs
            for atom in molecule.structure[start]:
                # only traverse through chain, ignore terminal atoms
                if not get_element(atom) in ["H", "F", "Cl", "Br", "I"]:
                    # no jumping forth and back between current and last visited node
                    if atom != ancestors[start]:
                        self._cycle_detection(molecule, atom, start, status, ancestors)
            # node done, mark node done
            status[start] = "KNOWN"


    # delete peptide bonds from lsit of rotatable bonds (torsions)
    def find_peptidebonds(self, molecule: Structure):
        peptidebonds = []

        # loop over every rotatable bond
        for bond in molecule.torsions:
            atom1 = bond[0]
            atom2 = bond[1]

            # check for characteristic C=O + C-NHR bonding environment
            # case first atom in list C, second atom in list N
            #
            # GEHT DAS ELEGANTER?? UND GENAUER, ZB MIT PRÜFUNG VON ANTEILIGER MERHFACHBINDUNG??
            element1 = get_element(atom1)
            element2 = get_element(atom2)
            # first check: C and N included in bond?
            if element1 == "C" and element2 == "N":
                # second check: C with bonding partner O?
                for atom3 in molecule.structure[atom1]:
                    element3 = get_element(atom3)
                    if element3 == "O":
                        # third check: double bond between C and O?
                        distance = np.sqrt((molecule.coords[atom1][0] - molecule.coords[atom3][0]) ** 2 +
                                           (molecule.coords[atom1][1] - molecule.coords[atom3][1]) ** 2 +
                                           (molecule.coords[atom1][2] - molecule.coords[atom3][2]) ** 2)

                        double_bond = covalence_radii[element1][1] + covalence_radii[element3][1] + 0.08

                        if distance <= double_bond:
                            # found peptide bond, remove bond from list of rotatable bonds
                            peptidebonds.append(bond)
                            break

            # check for characteristic C=O + C-NHR bonding environment
            # first atom in list C, second atom in list N
            elif element1 == "N" and element2 == "C":
                # first check: C and N included in bond?
                for atom3 in molecule.structure[atom2]:
                    element3 = get_element(atom3)
                    # second check: C with bonding partner O?
                    if element3 == "O":
                        # third check: double bond between C and O?
                        distance = np.sqrt((molecule.coords[atom2][0] - molecule.coords[atom3][0]) ** 2 +
                                           (molecule.coords[atom2][1] - molecule.coords[atom3][1]) ** 2 +
                                           (molecule.coords[atom2][2] - molecule.coords[atom3][2]) ** 2)

                        double_bond = covalence_radii[element2][1] + covalence_radii[element3][1] + 0.08

                        # found peptide bond, remove bond from list of rotatable bonds
                        if distance <= double_bond:
                            peptidebonds.append(bond)
                            break
        if peptidebonds:
            while True:
                rotatable = input("Consider peptide bonds as rotatable (1) or not (2)?: ")
                if rotatable == "1":
                    break
                elif rotatable == "2":
                    for bond in peptidebonds:
                        molecule.torsions.remove(bond)
                    break
                else:
                    print("Invalid input.")


    # delete peptide bonds from lsit of rotatable bonds (torsions)
    def __find_peptidebonds_old(self, molecule: Structure):
        iter = 0

        # loop over every rotatable bond
        while iter < len(molecule.torsions):
            atom1 = molecule.torsions[iter][0]
            atom2 = molecule.torsions[iter][1]

            # check for characteristic C=O + C-NHR bonding environment
            # case first atom in list C, second atom in list N
            #
            # GEHT DAS ELEGANTER?? UND GENAUER, ZB MIT PRÜFUNG VON ANTEILIGER MERHFACHBINDUNG??
            element1 = get_element(atom1)
            element2 = get_element(atom2)
            # first check: C and N included in bond?
            if element1 == "C" and element2 == "N":
                # second check: C with bonding partner O?
                for atom3 in molecule.structure[atom1]:
                    element3 = get_element(atom3)
                    if element3 == "O":
                        # third check: double bond between C and O?
                        distance = np.sqrt((molecule.coords[atom1][0] - molecule.coords[atom3][0]) ** 2 +
                                           (molecule.coords[atom1][1] - molecule.coords[atom3][1]) ** 2 +
                                           (molecule.coords[atom1][2] - molecule.coords[atom3][2]) ** 2)

                        double_bond = covalence_radii[element1][1] + covalence_radii[element3][1] + 0.08

                        if distance <= double_bond:
                            # found peptide bond, remove bond from list of rotatable bonds
                            molecule.torsions.remove(molecule.torsions[iter])
                            break

            # check for characteristic C=O + C-NHR bonding environment
            # first atom in list C, second atom in list N
            elif element1 == "N" and element2 == "C":
                # first check: C and N included in bond?
                for atom3 in molecule.structure[atom2]:
                    element3 = get_element(atom3)
                    # second check: C with bonding partner O?
                    if element3 == "O":
                        # third check: double bond between C and O?
                        distance = np.sqrt((molecule.coords[atom2][0] - molecule.coords[atom3][0]) ** 2 +
                                           (molecule.coords[atom2][1] - molecule.coords[atom3][1]) ** 2 +
                                           (molecule.coords[atom2][2] - molecule.coords[atom3][2]) ** 2)

                        double_bond = covalence_radii[element2][1] + covalence_radii[element3][1] + 0.08

                        # found peptide bond, remove bond from list of rotatable bonds
                        if distance <= double_bond:
                            molecule.torsions.remove(molecule.torsions[iter])
                            break
            iter += 1


    # control function for selection menu with thread for manual selection and graphical output
    def _selection_menu(self, molecule: Structure) -> bool:
        selection = [False]
        # if there rotatable bonds in the molecule start selection menu
        if molecule.torsions:
            # start manual selection
            selection_menu = threading.Thread(target=self._select_torsions, args=[molecule, selection], daemon=True)
            selection_menu.start()
            # start graphical output
            self.show_structure_colors(molecule, selection_menu)
        return selection[0]


    # graphical output of the input molecule structure with numbering of the atoms
    def show_structure(self, molecule: Structure, selection_thread: threading.Thread = None):
        # lists for coordinates, bonds and labels of chain atoms
        central_xcoords = list()
        central_ycoords = list()
        central_zcoords = list()
        central_bonds = list()
        central_labels = list()

        # lists for coordinates, bonds and labels of terminal atoms
        terminal_xcoords = list()
        terminal_ycoords = list()
        terminal_zcoords = list()
        terminal_bonds = list()
        terminal_labels = list()

        # loop over all atoms
        for atom1 in molecule.structure:
            # get coordinates of first atom of a bond
            xcoord1 = molecule.coords[atom1][0]
            ycoord1 = molecule.coords[atom1][1]
            zcoord1 = molecule.coords[atom1][2]

            # if atom is terminal raise flag
            terminal1 = get_element(atom1) in ["H", "F", "Cl", "Br", "I"]

            # if flag is raised store coordinates and label of atom in according lists of terminal atoms
            if terminal1:
                terminal_xcoords.append(xcoord1)
                terminal_ycoords.append(ycoord1)
                terminal_zcoords.append(zcoord1)

                terminal_labels.append(atom1)
            # if flag is raised store coordinates and label of atom in according lists of chain atoms
            else:
                central_xcoords.append(xcoord1)
                central_ycoords.append(ycoord1)
                central_zcoords.append(zcoord1)

                central_labels.append(atom1)

            # loop over every bond partner of atom
            for atom2 in molecule.structure[atom1]:
                # get coordinates of bond partner
                xcoord2 = molecule.coords[atom2][0]
                ycoord2 = molecule.coords[atom2][1]
                zcoord2 = molecule.coords[atom2][2]

                # if bond partner is terminal raise flag
                terminal2 = get_element(atom2) in ["H", "F", "Cl", "Br", "I"]

                # create bond for atom and bond partner in according list
                new_bond = [[xcoord1, xcoord2], [ycoord1, ycoord2], [zcoord1, zcoord2]]

                # if atom or/and bond partner is terminal store bond in list of terminal bonds
                if terminal1 or terminal2:
                    terminal_bonds.append(new_bond)
                # if both atom and bond partner are not terminal store bond in list of central bonds
                else:
                    central_bonds.append(new_bond)

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

        # plot central bonds (black)
        for bond in central_bonds:
            ax.plot([bond[x][atom1], bond[x][atom2]], [bond[y][atom1], bond[y][atom2]],
                    [bond[z][atom1], bond[z][atom2]], color="black")

        # plot terminal bonds (light gray)
        for bond in terminal_bonds:
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


    # graphical output of the input molecule structure with numbering of the atoms and coloring of bonds
    def show_structure_colors(self, molecule: Structure, selection_thread: threading.Thread = None):
        # lists for coordinates, bonds and labels of chain atoms
        central_xcoords = list()
        central_ycoords = list()
        central_zcoords = list()
        central_bonds1 = list()
        central_bonds2 = list()
        central_bonds3 = list()
        central_labels = list()

        # lists for coordinates, bonds and labels of terminal atoms
        terminal_xcoords = list()
        terminal_ycoords = list()
        terminal_zcoords = list()
        terminal_bonds = list()
        terminal_labels = list()

        # loop over all atoms
        for atom1 in molecule.structure:
            # get coordinates of first atom of a bond
            xcoord1 = molecule.coords[atom1][0]
            ycoord1 = molecule.coords[atom1][1]
            zcoord1 = molecule.coords[atom1][2]

            # if atom is terminal raise flag
            terminal1 = get_element(atom1) in ["H", "F", "Cl", "Br", "I"]

            # if flag is raised store coordinates and label of atom in according lists of terminal atoms
            if terminal1:
                # store coordinates
                terminal_xcoords.append(xcoord1)
                terminal_ycoords.append(ycoord1)
                terminal_zcoords.append(zcoord1)

                # store label
                terminal_labels.append(atom1)
            # if flag is raised save coordinates and label of atom in according lists of chain atoms
            else:
                # store coordinates
                central_xcoords.append(xcoord1)
                central_ycoords.append(ycoord1)
                central_zcoords.append(zcoord1)

                # store label
                central_labels.append(atom1)

            # loop over every bond partner of atom
            for atom2 in molecule.structure[atom1]:
                # get coordinates of bond partner
                xcoord2 = molecule.coords[atom2][0]
                ycoord2 = molecule.coords[atom2][1]
                zcoord2 = molecule.coords[atom2][2]

                # if bond partner is terminal raise flag
                terminal2 = get_element(atom2) in ["H", "F", "Cl", "Br", "I"]

                # create bond for atom and bond partner in according list
                new_bond = [[xcoord1, xcoord2], [ycoord1, ycoord2], [zcoord1, zcoord2]]

                # if atom or/and bond partner is terminal store bond in list of terminal bonds
                if terminal1 or terminal2:
                    terminal_bonds.append(new_bond)

                # if bond does not include terminal atoms (central bond) and is not rotatable store in list
                # for concerning coloring in plot
                elif not ([atom1, atom2] in molecule.torsions or [atom2, atom1] in molecule.torsions):
                    central_bonds3.append(new_bond)

                # if bond does not include terminal atoms (central bond) and is rotatable
                else:
                    # count non-terminal neighbors of first atom in bond
                    neighbors1 = 0
                    for neighbor in molecule.structure[atom1]:
                        if not get_element(neighbor) in ["H", "F", "Cl", "Br", "I"]:
                            neighbors1 += 1

                    # if first atom of bond has at least two non-terminal neighbors bond could be central, if not
                    # it can already be declared as bond to a terminal group and stored in concerning list
                    if neighbors1 >= 2:
                        # count non-terminal neighbors of second in atom in bond
                        neighbors2 = 0
                        for neighbor in molecule.structure[atom2]:
                            if not get_element(neighbor) in ["H", "F", "Cl", "Br", "I"]:
                                neighbors2 += 1

                        # if second atom of bond has at least two non-terminal neighbors bond is central and stored
                        # in list concerning coloring in plot, if not store it in list with other color
                        if neighbors2 >= 2:
                            central_bonds2.append(new_bond)
                        else:
                            central_bonds1.append(new_bond)
                    else:
                        central_bonds1.append(new_bond)

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
        for bond in central_bonds1:
            ax.plot([bond[x][atom1], bond[x][atom2]], [bond[y][atom1], bond[y][atom2]],
                    [bond[z][atom1], bond[z][atom2]], color="blue")

        # plot rotatable, central bonds in red
        for bond in central_bonds2:
            ax.plot([bond[x][atom1], bond[x][atom2]], [bond[y][atom1], bond[y][atom2]],
                    [bond[z][atom1], bond[z][atom2]], color="red")

        # plot non-rotatable bonds in black
        for bond in central_bonds3:
            ax.plot([bond[x][atom1], bond[x][atom2]], [bond[y][atom1], bond[y][atom2]],
                    [bond[z][atom1], bond[z][atom2]], color="black")

        # plot terminal bonds in light gray
        for bond in terminal_bonds:
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


    # menu for user to select bonds to be rotated around and rotation angle increment
    def _select_torsions(self, molecule: Structure, selection: list):
        # check1 controls if menu is repeated or closed after selection
        check1 = False
        while not check1:
            #selection = False
            # default values for rotation angle increment and list of bonds to be rotated around
            self.angle_increments = 0
            selected_torsions = list()
            # show list with rotatable bonds and selection options
            print("Rotatable bonds:")
            for number, torsion in enumerate(molecule.torsions):
                print("Nr." + str(number) + ": " + torsion[0] + "-" + torsion[1])
            print("Vary torsion angles of all bonds, which means blue+red (1)")
            print("or vary only torsion angles of central bonds, which means blue (2)")
            print("or select bonds manually (3)")

            # Bei sinnvoller Eingabe (1,2 oder 3) gilt Block als erledigt (check2=True), bei nicht sinnvoller Eingabe
            # Wiederholung des Blocks (check2=False)
            # check2 controls if mode selection is repeated or not
            check2 = False
            while not check2:
                mode_input = input("Type in number and press Enter: ")

                # mode: all bonds (1)
                if mode_input == "1":
                    print()
                    # mode selection finished
                    check2 = True
                    # selection of rotation angle increment
                    self._get_angle_increment()
                    # confirm choices and finish menu or not
                    check1 = self._get_conformation(molecule)
                    selection[0] = True

                # mode: just central bonds (2)
                elif mode_input == "2":
                    print()
                    # mode selection finished
                    check2 = True

                    # loop over all rotatable bonds, check if bond is central
                    for torsion in molecule.torsions:
                        atom1 = torsion[0]
                        atom2 = torsion[1]

                        # count chain/non-terminal atoms which first atom of bond is bound to
                        chain_neighbors1 = 0
                        for neighbor in molecule.structure[atom1]:
                            if not get_element(neighbor) in ["H", "F", "Cl", "Br", "I"]:
                                chain_neighbors1 += 1
                        # if first atom is bound to at least two chain/non-terminal atoms bond could be central
                        # count chain/non-terminal atoms which second atom of bond is bound to
                        if chain_neighbors1 >= 2:
                            chain_neighbors2 = 0
                            for neighbor in molecule.structure[atom2]:
                                if not get_element(neighbor) in ["H", "F", "Cl", "Br", "I"]:
                                    chain_neighbors2 += 1
                            # if also second atom is bound to at least two chain/non-terminal atoms bond is central
                            # and is therefore included in list of bonds to be rotated
                            if chain_neighbors2 >= 2:
                                selected_torsions.append(torsion)

                    # selection of rotation angle increment
                    self._get_angle_increment()
                    # confirm choices and finish menu or not
                    check1 = self._get_conformation(molecule, selected_torsions)
                    selection[0] = True

                # mode: manual bond selection (3)
                elif mode_input == "3":
                    print()
                    # mode selection finished
                    check2 = True
                    print("Type in the numbers of the bonds which should be rotated. Finish by typing in -1.")

                    # bonds which can be chosen (rotatable bonds), once a bond is selected it is removed from list
                    # to handle wrong input (double selection, non-existing bond, selection of no bond)
                    user_selection = [i for i in range(len(molecule.torsions))]

                    # ask for input as long as there are bonds left which can be chosen
                    while True and user_selection:
                        bond_index = input("Bond: ")
                        # input -1: selection finished
                        if bond_index == "-1":
                            # no bond has been choosen
                            if len(user_selection) == len(molecule.torsions):
                                while True:
                                    confirmation = input("No bond number has been given, there will no conformers be "
                                                         "generated. Finish (1) or continue selecting (2)?: ")
                                    # end of menu
                                    if confirmation == "1":
                                        selection[0] = False
                                        return
                                    # continue selection
                                    elif confirmation == "2":
                                        break
                                    # wrong input, repeat asking for confirmation
                                    else:
                                        print("Invalid input.")
                            # at least one bond has been chosen, end bond selection
                            else:
                                break
                        # invalid input (letter, float, negative number), repeat asking for input
                        elif not bond_index.isdecimal():
                            print("Invalid input.")
                            continue
                        # invalid bond number, repeat asking for input
                        elif int(bond_index) >= len(molecule.torsions):
                            print("Bond number doesn't exist (see above).")
                            continue
                        # valid input/bond number
                        else:
                            bond_index = int(bond_index)
                            # check if bond has already been selected, if so repeat asking for input
                            try:
                                user_selection.remove(bond_index)
                            except ValueError:
                                print("Bond has already been selected.")
                                continue
                            # if bond has not been selected yet store in list of selected bonds
                            selected_torsions.append(molecule.torsions[bond_index])
                    # selection of rotation angle increment
                    self._get_angle_increment()
                    # confirm choices and finish menu or not
                    check1 = self._get_conformation(molecule)
                    selection[0] = True
                # invalid input in mode selection, repeat asking for input
                else:
                    print("Invalid input.")
                    continue


    # user chooses angle increment for variation of dihedral angles in generation of conformers structures
    def _get_angle_increment(self):
        while True:
            # get input for rotation angle inrement
            increment = input("Type in and confirm the angle increment by which the torsion angles will be varied"
                              " (between 0 and 360 degrees)\nor press Enter to accept default (120 degrees): ")
            # check if input is a number
            if increment.isdigit():
                # increments higher than 360 degrees or lower than 0 degrees are turned into equivalent
                # angle increments between 0 and 360 degrees
                increment = int(increment) % 360
                # 0 as increment is invalid to prevent division by zero in some other function
                if increment != 0:
                    self.angle_increment = increment
                    return
                else:
                    print("Invalid input.")
            elif increment == "":
                return
            else:
                print("Invalid input.")


    # Bestätigung der Auswahl durch User
    #
    # WAS WENN WINKEL 0 ANGEGEBEN WIRD??
    def _get_conformation(self, molecule: Structure, selected_torsions: list = None) -> bool:
        # je nach vorherige Auswahl der Bindungen speichern der Anzahl der zu rotierenden Bindungen für
        # Berechnung der resultierenden Konformerenanzahl
        if selected_torsions:
            torsions = len(selected_torsions)
        else:
            torsions = len(molecule.torsions)

        while True:
            # Berechnung und Ausgabe von Anzahl der zu generierenden Konformere, Abfrage ob damit Fortsetzen (1)
            # oder neue Auswahl (2)
            amount = int(np.power((360 / self.angle_increment), torsions))
            print()
            confirmation = input("Up to " + str(amount) + " conformers will be generated. Continue (1) or start"
                                                          " new selection (2)?: ")

            # Eingabe 1: Fortsetzen, in selected_torsions zwischengespeicherte Bindungen werden in torsions
            # übernommen dieser Block (check4=True) und der übergeordnete Block zur Modus-Abfrage (check1=True)
            # gelten als erledigt
            if confirmation == "1":
                print()
                if selected_torsions:
                    molecule.torsions = selected_torsions
                return True

            # Eingabe 2: Neue Auswahl, dieser Block gilt als erledigt, der übergeordnete
            # Block zur Modus-Abfrage wird wiederholt (check1=False, nächster Schleifendurchlauf)
            elif confirmation == "2":
                print()
                return False


    def generate_conformers(self, molecule: Structure, constraint_opt: bool = False, cluster: bool = False):
        if not molecule.coords and not molecule.structure:
            print("CONFORMERGENERATOR MODULE: Molecule coordinates and structural information are missing.")
            return
        if not molecule.coords:
            print("CONFORMERGENERATOR MODULE: Molecule coordinates are missing.")
            return
        if not molecule.structure:
            print("CONFORMERGENERATOR MODULE: Molecule structural information is missing.")
            return

        self.find_cycles(molecule)
        self.find_peptidebonds(molecule)
        selection = self._selection_menu(molecule)

        if selection:
            if cluster:
                #HARDCODING VERMEIDEN
                self.angles = [0, 120, 240]
                self._combinations(molecule)
            else:
                if self.angle_increment:
                    self._generation_setup(molecule)

                    if constraint_opt:
                        conformers = self._combinations_constraints(molecule)
                    else:
                        conformers = self._combinations(molecule)

                    if conformers:
                        print(str(conformers) + " conformers have been generated.")
                    else:
                        print("No conformers could be generated due to clashes in every calculated structure.")

        """
        if cluster:
            #HARDCODING VERMEIDEN
            self.angles = [0, 120, 240]
            self.combinations(molecule)
        else:
            if self.angle_increment:
                self.generation_setup(molecule)

                if constraint_opt:
                    conformers = self.combinations_constraints(molecule)
                else:
                    conformers = self.combinations(molecule)

                if conformers:
                    print(str(conformers) + " conformers have been generated.")
                else:
                    print("No conformers could be generated due to clashes in every calculated structure.")
        """


    # zählen der jeweils zu rotierenden Atome auf beiden Seiten der Bindung, Hinzufügen zu torsion_atoms von der
    # Seite welche weniger Atome enthält
    def _generation_setup(self, molecule: Structure):
        self.angles = [n * self.angle_increment for n in range(int(np.round(360 / self.angle_increment)))]

        # Durchführung für jede gefundene und gefilterte bzw. ausgewählte Rotationsbindung
        for bond_number, bond in enumerate(molecule.torsions):
            self.torsion_atoms.append([])
            self.torsion_atoms.append([])

            # zählen der linken Seite (entspricht erstgenanntem Atom der Bindung in torsions)
            status = {atom: "UNKNOWN" for atom in molecule.structure}
            left = self._torsion_atom_counter(molecule, bond[0], bond[1], status, -1, bond_number)

            # zählen der linken Seite (entspricht letztgenanntem Atom der Bindung in torsions)
            status = {atom: "UNKNOWN" for atom in molecule.structure}
            right = self._torsion_atom_counter(molecule, bond[1], bond[0], status, -1, bond_number + 1)

            # wenn links weniger oder genauso viele Atome wie rechts wird links gedreht, d.h. Atome der linken
            # Seite zu torsions_atoms hinzufügen
            if left <= right:
                self.torsion_atoms[bond_number].remove(bond[0])
                self.torsion_atoms.remove(self.torsion_atoms[bond_number + 1])
            # wenn rechts weniger  Atome wie links wird links gedreht, d.h. Atome der rechten Seite zu
            # torsions_atoms hinzufügen
            else:
                self.torsion_atoms[bond_number + 1].remove(bond[1])
                self.torsion_atoms.remove(self.torsion_atoms[bond_number])


    # Hilfsfunktion, wird von rotation_setup aufgerufen und zählt zu rotierende Atome ausgehend vom Atom einer
    # Bindung bei Rotation um diese
    #
    # basiert auf Depth-First-Search (DFS) auf ungerichtetem Graphen mit Färbemethode (analog zu find_cycles)
    # UNKNOWN: Knoten noch nicht besucht
    # SEEN: Knoten besucht
    def _torsion_atom_counter(self, molecule: Structure, atom, last_atom, status, counter, bond):
        # Atom bereits betrachtet
        if status[atom] == "SEEN":
            return counter
        # Atom ist Terminalatom, d.h. Ende einer Kette
        elif get_element(atom) in ["H", "F", "Cl", "Br", "I"]:
            # Atom als betrachtet markieren
            status[atom] = "SEEN"
            # Atom zu torsion_atoms hinzufügen
            self.torsion_atoms[bond].append(atom)
            return counter + 1
        # Atom noch nicht betrachtet und nicht Ende einer Kette
        else:
            # Atom als betrachtet markieren
            status[atom] = "SEEN"
            # Atom zu torsion_atoms hinzufügen
            self.torsion_atoms[bond].append(atom)
            # Zähler aktualisieren
            counter += 1
            # rekursiver Aufruf der Funktion für alle Bindungspartner von Atom, solange diesem nicht dem
            # Atom des vorherigen Funktionsaurufs entsprechen
            for neighbor in molecule.structure[atom]:
                if neighbor != last_atom:
                    counter = self._torsion_atom_counter(molecule, neighbor, atom, status, counter, bond)
        # Anzahl der zu rotierenden Atome der entsprechenden Seite der Bindung zurückgeben
        return counter


    # calculates all possible conformer structures and generates an output file for every conformer structure without
    # internal clash
    # the output file is an input file for an geometry optimization with ORCA
    def _combinations(self, molecule: Structure, new_coords: dict = None, index: int = 0, counter: int=0) -> int:
        # base case, new torsion angle for every angle has been calculated
        if index == len(molecule.torsions):
            # sofern keine strukturinternen Clashes hinzufügen zur Liste erfolgreich erzeugter Konformerstrukturen
            # check new structure for internal clashes
            if not molecule._clashes(new_coords):
                self.output_opt(new_coords, counter)
                return counter+1
            else:
                return counter

        # es wurden noch nicht alle Torsionswinkel berechnet, Bindung index in torsions ist an der Reihe
        else:
            if not new_coords:
                new_coords = molecule.coords

            # Punkte initialisieren, welche die Rotationsachse definieren
            vec1 = np.array(new_coords[molecule.torsions[index][0]])
            vec2 = np.array(new_coords[molecule.torsions[index][1]])

            # jeden möglichen Torsionswinkel für Bindung durchgehen
            for angle in self.angles:
                # Python übergibt quasi grundsätzlich per Referenz, mit Kopie arbeiten damit neue Winkel für
                # jede Bindung ausgehend von unverändertem Torsionswinkel in dieser Bindung berechnet werden
                #
                # ODER AUCH SCHON VOR FOR-SCHLEIFE NÖTIG??
                new_coords_copy = new_coords.copy()

                # Rotation aller zu Bindung gehöriger Atome in torsions_atoms um betrachtete Rotationsachse
                #
                # VERBESSERUNG WOMÖGLICH WENN ALLE PUNKTE AUF EINMAL BERECHNET WERDEN (MIT DIESEM "TILE" UND SO)
                for atom in self.torsion_atoms[index]:
                    coords = np.array(new_coords[atom])
                    coords = rotation(vec2, vec1, angle, coords)
                    new_coords_copy[atom] = list(coords)

                # rekursiver Aufruf für nächste Bindung
                counter = self._combinations(molecule, new_coords_copy, index+1, counter)
            return counter


    # calculates all possible conformer structures and generates an output file for every conformer structure without
    # internal clash
    # the output file is an input file for an constraint geometry optimization with ORCA
    def _combinations_constraints(self, molecule: Structure, new_coords: dict = None,
                                 index: int = 0, counter: int = 0, constraints: list = None) -> int:
        # base case, Torsionswinkel für jede Bindung wurde berechnet
        if index == len(molecule.torsions):
            # sofern keine strukturinternen Clashes hinzufügen zur Liste erfolgreich erzeugter Konformerstrukturen
            if not molecule._clashes(new_coords):
                #self.conformers.append(new_coords)
                self.output_constraint_opt(new_coords, constraints, counter)
                return counter + 1
            else:
                return counter

        # es wurden noch nicht alle Torsionswinkel berechnet, Bindung index in torsions ist an der Reihe
        else:
            if not new_coords:
                new_coords = molecule.coords
                constraints = []

            # Punkte initialisieren, welche die Rotationsachse definieren
            vec1 = np.array(new_coords[molecule.torsions[index][0]])
            vec2 = np.array(new_coords[molecule.torsions[index][1]])

            # jeden möglichen Torsionswinkel für Bindung durchgehen
            for angle in self.angles:
                # Python übergibt quasi grundsätzlich per Referenz, mit Kopie arbeiten damit neue Winkel für
                # jede Bindung ausgehend von unverändertem Torsionswinkel in dieser Bindung berechnet werden
                #
                # ODER AUCH SCHON VOR FOR-SCHLEIFE NÖTIG??
                new_coords_copy = new_coords.copy()
                # Block für beschränkte Optimierungen
                curr_constraints = constraints.copy()
                if angle != 0:
                    neighbor_index = 0
                    axis_atom1 = molecule.torsions[index][0]
                    axis_atom2 = molecule.torsions[index][1]
                    neighbor1 = molecule.structure[axis_atom1][neighbor_index]
                    while neighbor1 == axis_atom1 or neighbor1 == axis_atom2:
                        neighbor_index += 1
                        neighbor1 = molecule.structure[axis_atom1][neighbor_index]
                    neighbor_index = 0
                    neighbor2 = molecule.structure[axis_atom2][neighbor_index]
                    while neighbor2 == axis_atom1 or neighbor2 == axis_atom2:
                        neighbor_index += 1
                        neighbor2 = molecule.structure[axis_atom2][neighbor_index]
                    curr_constraints.append([get_number(neighbor1), get_number(axis_atom1),
                                             get_number(axis_atom2), get_number(neighbor2)])

                # Rotation aller zu Bindung gehöriger Atome in torsions_atoms um betrachtete Rotationsachse
                #
                # VERBESSERUNG WOMÖGLICH WENN ALLE PUNKTE AUF EINMAL BERECHNET WERDEN (MIT DIESEM "TILE" UND SO)
                for atom in self.torsion_atoms[index]:
                    coords = np.array(new_coords[atom])
                    coords = rotation(vec2, vec1, angle, coords)
                    new_coords_copy[atom] = list(coords)

                # rekursiver Aufruf für nächste Bindung
                counter = self._combinations_constraints(molecule, new_coords_copy, index + 1,
                                                        counter, curr_constraints)
            return counter


    def output_coords(self, counter: int):
        # Fall: es wurde mindestens ein Konformer erzeugt
        conformers = len(list(self.conformers))
        if conformers:
            if not os.path.exists("../Output"):
                os.makedirs("../Output")
            for counter, conformer in enumerate(self.conformers):
                # schreiben der Koordinaten in Outputdatei
                #
                # ALLGEMEINE FUNKTION IN HELPER DARAUS MACHEN??
                with open("../Output/conformer" + str(counter) + ".xyz", "w") as outfile:
                    # outfile.write(str(len(new_coords)))
                    # outfile.write("\n\n")
                    for atom in conformer:
                        outfile.write(get_element(atom))
                        outfile.write("\t")
                        print("{:18.15f} \t {:18.15f} \t {:18.15f}".format(conformer[atom][0], conformer[atom][1],
                                                                           conformer[atom][2]), file=outfile)
            print(str(conformers) + " conformers have been generated.")

        # Fall: es wurde wegen interner Clashes bei jeder generierten Struktur kein Konformer erzeugt
        else:
            print("No conformers could be generated because of clashes in every calculated structure.")


    def output_opt(self, coords: dict, counter: int):
        if not os.path.exists("../Output"):
            os.makedirs("../Output")
        with open("../Output/conformer" + str(counter), "w") as outfile:
            outfile.write("! BP86 SVP D4 TightSCF TightOpt NumFreq PAL4\n\n")
            outfile.write("*xyz 0 1\n")
            for atom in coords:
                outfile.write(get_element(atom))
                outfile.write("\t")
                print("{:18.15f} \t {:18.15f} \t {:18.15f}".format(coords[atom][0], coords[atom][1],
                                                                   coords[atom][2]), file=outfile)
            outfile.write("*")


    def output_constraint_opt(self, coords: dict, constraints: list, counter: int):
        if not os.path.exists("../Output"):
            os.makedirs("../Output")
        with open("../Output/conformer" + str(counter), "w") as outfile:
            outfile.write("! BP86 SVP D4 TightSCF TightOpt NumFreq PAL4\n\n")
            if constraints:
                outfile.write("%geom Constraints \n")
                for constraint in constraints:
                    print("{{D {:d} {:d} {:d} {:d} C}}".format(constraint[0], constraint[1],
                                                                   constraint[2], constraint[3]), file=outfile)
                outfile.write("end\nend\n\n")
            outfile.write("*xyz 0 1\n")
            for atom in coords:
                outfile.write(get_element(atom))
                outfile.write("\t")
                print("{:18.15f} \t {:18.15f} \t {:18.15f}".format(coords[atom][0], coords[atom][1],
                                                                   coords[atom][2]), file=outfile)
            outfile.write("*")
