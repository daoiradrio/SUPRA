from utils.helper import is_hb_don, is_hb_acc, rotation, get_element, get_number
from SUPRACluster.clusterstructure import ClusterStructure
from SUPRAConformer.conformergenerator import ConformerGenerator
import numpy as np
import os
from copy import deepcopy



# TODO
# - WOZU WIRD MONOMER_STRUCTURE GESPEICHERT?? WIRD NICHT OHNEHIN STETS EINE CLUSTERSTRUCTURE MIT ÜBERGEBEN
#   WODURCH WEITERES SPEICHERN DER KONNEKTIVITÄT IM GENERATOR UNNÖTIG IST??
#   --> WENN DAS ANGEPASST WIRD MUSS AUCH HANDLING DER Z-MATRIZEN ÜBERARBEITET WERDEN



class ClusterGenerator:

    #SO KEINE FLEXIBILTÄT BZGL. CLUSTERN AUS VERSCHIEDENARTIGEN MONOMEREN (S.O. TO DO)
    #HIER MUSS ANGEPASST WERDEN WENN VERSCHIEDENARTIGE MONOMERE VERBAUT WERDEN (self.monomer_structure)
    def __init__(self, monomer: ClusterStructure):
        self.hb_len = 1.1
        self.zmatrices = dict()
        self.container = list()
        self.monomer = monomer
        self.counter = 0
        self.max_cluster_size = 3 # VOM USER WÄHLEN LASSEN
    

    def add_monomer(self, cluster1: ClusterStructure, atom1: str, cluster2: ClusterStructure, atom2: str):
        if is_hb_don(atom1):
            if not is_hb_acc(atom2):
                print("ERROR")
                return
            docking_cluster = cluster1
            docking_atom = atom1
            dock_cluster = cluster2
            dock_atom = atom2
        elif is_hb_don(atom2):
            if not is_hb_acc(atom1):
                print("ERROR")
                return
            docking_cluster = cluster2
            docking_atom = atom2
            dock_cluster = cluster1
            dock_atom = atom1
        label_shift = len(dock_cluster.coords.keys())
        pos = dock_cluster.get_acc_vec(dock_atom)
        coords_shift = pos - docking_cluster.coords[docking_atom]
        v1 = dock_cluster.coords[dock_atom] - pos
        v1 = v1 / np.linalg.norm(v1)
        v2 = docking_cluster.get_don_vec(docking_atom) + coords_shift - pos
        v2 = v2 / np.linalg.norm(v2)
        angle = np.arccos(np.dot(v1, v2))
        norm_vec = np.cross(v2, v1) + pos
        for atom, coords in docking_cluster.coords.items():
            new_label = f"{get_element(atom)}{get_number(atom) + label_shift}"
            new_coords = coords + coords_shift
            new_coords = rotation(new_coords, pos, norm_vec, angle)
            dock_cluster.coords[new_label] = new_coords


    def new_set_zmatrix(self, hb_atom: str):
        new_zmatrix = []
        self.zmatrices[hb_atom] = new_zmatrix
        if self.monomer.get_element(hb_atom) == "H":
            atom1 = hb_atom
            atom2 = self.monomer.bond_partners[hb_atom][0]
            pos1 = self.monomer.coords[atom1]
            pos2 = self.monomer.coords[atom2]
            distance = np.linalg.norm(pos2 - pos1)
            new_zmatrix.append([atom1])
            new_zmatrix.append([atom2, distance])
            atoms = [atom for atom in self.monomer.coords.keys() if atom != hb_atom and atom != atom2]
        elif self.monomer.get_element(hb_atom) in ["O", "N"]:
            atom2 = hb_atom
            pos1 = self.monomer.get_don_vec(hb_atom)
            pos2 = self.monomer.coords[atom2]
            new_zmatrix.append([atom2, 1.0])
            atoms = [atom for atom in self.monomer.coords.keys() if atom != hb_atom]
        atom3 = atoms[0]
        pos3 = self.monomer.coords[atom3]
        distance = np.linalg.norm(pos3 - pos2)
        angle = self.get_angle(pos1, pos2, pos3)
        new_zmatrix.append([atom3, distance, angle])
        for atom4 in atoms[1:]:
            pos4 = self.monomer.coords[atom4]
            distance = np.linalg.norm(pos4 - pos3)
            angle = self.get_angle(pos2, pos3, pos4)
            dihedral = self.get_dihedral(pos1, pos2, pos3, pos4)


    def get_angle(self, p1: np.array, p2: np.array, p3: np.array):
        print(p1, p2, p3)
        v1 = p2 - p1 / np.linalg.norm(p2 - p1)
        v2 = p3 - p2 / np.linalg.norm(p3 - p2)
        cos = np.dot(v1, v2)
        sin = np.linalg.norm(np.cross(v1, v2))
        return np.arctan2(sin, cos)


    def get_dihedral(self, p1: np.array, p2: np.array, p3: np.array, p4: np.array) -> None:
        v1 = p2 - p1 / np.linalg.norm(p2 - p1)
        v2 = p3 - p2 / np.linalg.norm(p3 - p2)
        v3 = p4 - p3 / np.linalg.norm(p4 - p3)
        n1 = np.cross(v2, v1)
        n1 = n1 / np.linalg.norm(n1)
        n2 = np.cross(v3, v1)
        n2 = n2 / np.linalg.norm(n2)


    def old_get_angle_math_stack_exchange(self, p1: np.array, p2: np.array, p3: np.array) -> float:
        v1 = (p2 - p1) / np.linalg.norm(p2 - p1)
        v2 = (p3 - p2) / np.linalg.norm(p3 - p2)
        norm_a = np.linalg.norm(v1 - v2)
        norm_b = np.linalg.norm(v1 + v2)
        return np.arctan2(norm_a, norm_b)


    def old_get_angle_see_ref_on_laptop(self, p1: np.array, p2: np.array, p3: np.array) -> float:
        v1 = p2 - p1
        v2 = p3 - p2
        v3 = v1 + v2
        a = np.linalg.norm(v1)
        b = np.linalg.norm(v2)
        c = np.linalg.norm(v3)
        a = max(a, b)
        b = max(a, b)
        if b >= c:
            y = c - (a - b)
        elif c >= b:
            y = b - (a - c)
        angle = 2 * np.arctan(
                    np.sqrt(
                        ((a - b) + c) * y / ((a + (b + c)) * ((a - c) + b))
                    )
                )
        return angle


    def set_zmatrix(self, hb_atom: str):
        new_zmatrix = []
        self.zmatrices[hb_atom] = new_zmatrix
        # case acceptor: set first (hb_atom, just label) and second atom (label & distance)
        if self.monomer.get_element(hb_atom) == "H":
            new_zmatrix.append([hb_atom])
            atom2 = self.monomer.bond_partners[hb_atom][0]
            vec12 = self.monomer.coords[atom2] - self.monomer.coords[hb_atom]
            len_vec12 = np.linalg.norm(vec12)
            vec12 = vec12 / len_vec12
            new_zmatrix.append([atom2, len_vec12])
            atoms = [atom for atom in self.monomer.coords.keys() if atom != hb_atom and atom != atom2]
        # case donator: set second atom (hb_atom, label & distance), first atom (just label) is skipped
        elif self.monomer.get_element(hb_atom) in ["O", "N"]:
            atom2 = hb_atom
            vec12 = self.monomer.get_don_vec(hb_atom) - self.monomer.coords[hb_atom]
            vec12 = vec12 / np.linalg.norm(vec12)
            new_zmatrix.append([hb_atom, 1.0]) #DISTANZ SO LASSEN?? MUSS JE NACH ELEMENTE KOMBINATION SOWIESO ANGEPASST WERDEN...
            atoms = [atom for atom in self.monomer.coords.keys() if atom != hb_atom] 
        # safety option for hydrogen bond elements implemented in ClusterStructure but not in ClusterGenerator yet
        else:
            print("Element not implemented here yet")
            return
        # set third atom
        atom3 = atoms[0]
        vec23 = self.monomer.coords[atom3] - self.monomer.coords[atom2]
        len_vec23 = np.linalg.norm(vec23)
        vec23 = vec23 / len_vec23
        angle = np.arccos(np.dot(-vec12, vec23))
        new_zmatrix.append([atom3, len_vec23, angle])
        ####################################
        print(hb_atom)
        print(self.get_angle(self.monomer.coords[hb_atom], self.monomer.coords[atom2], self.monomer.coords[atom3]))
        print(angle)
        print()
        ####################################
        # set remaining n-3 atoms
        for atom4 in atoms[1:]:
            vec34 = self.monomer.coords[atom4] - self.monomer.coords[atom3]
            len_vec34 = np.linalg.norm(vec34)
            vec34 = vec34 / len_vec34
            angle = np.arccos(np.dot(-vec23, vec34))
            # Wikipedia simplified formula for dihedral
            # https://en.wikipedia.org/wiki/Dihedral_angle#:~:text=A%20dihedral%20angle%20is%20the,line%20as%20a%20common%20edge.
            # u1 = vec12, u2 = vec23, u3 = vec34
            u1u2 = np.cross(vec12, vec23)
            u1u2 = u1u2 / np.linalg.norm(u1u2)
            u2u3 = np.cross(vec23, vec34)
            u2u3 = u2u3 / np.linalg.norm(u2u3)
            len_u2 = len_vec23
            dihedral = np.arctan2(len_u2 * np.dot(vec12, u2u3), np.dot(u1u2, u2u3))
            new_zmatrix.append([atom4, len_vec34, angle, dihedral])

    
    # case acceptor atom: atom1 = neighbor of acceptor, atom0 = acceptor
    # case donator atom: atom1 = donator
    def old_set_zmatrix(self, monomer: ClusterStructure, atom1: str, atom0: str=None):
        new_zmatrix = list()
        atom0_pos = []

        # add first atom to z-matrix
        # case acceptor
        if atom0:
            self.zmatrices[atom0] = new_zmatrix
            atoms_list = [atom for atom in monomer.coords if atom != atom0 and atom != atom1]
            atom0_pos = monomer.coords[atom0]
            new_zmatrix.append([atom0])
        # case donator
        else:
            self.zmatrices[atom1] = new_zmatrix
            atoms_list = [atom for atom in monomer.coords if atom != atom1]
            atom0_pos = monomer.coords[atom1] + 1.1 * monomer.hb_don_vec[get_number(atom1)]

        # add second atom to z-matrix
        atom1_pos = monomer.coords[atom1]
        vec10 = atom0_pos - atom1_pos
        distance10 = np.sqrt(np.dot(vec10, vec10))
        vec10 = vec10 / distance10
        new_zmatrix.append([atom1, distance10])

        #ERST ÜBERPRÜFEN OB ÜBERHAUPT EIN DRITTES ATOM VORHANDEN IST
        # add third atom to z-matrix (first atom from atoms_list)
        atom2 = atoms_list[0]
        atom2_pos = monomer.coords[atom2]
        vec12 = atom2_pos - atom1_pos
        distance = np.sqrt(np.dot(vec12, vec12))
        vec12 = vec12 / distance
        angle = np.arccos(np.dot(vec10, vec12))
        new_zmatrix.append([atom2, distance, angle])

        # if atoms_list contains 2 or more atoms add them to z-matrix
        for atom3 in atoms_list[1:]:
            # necessary initialization
            atom3_pos = monomer.coords[atom3]

            # distance between atom n and atom n-1
            vec23 = atom3_pos - atom2_pos
            distance = np.sqrt(np.dot(vec23, vec23))
            vec23 = vec23 / distance

            # angle between vector(atom1->atom2) and vector(atom2->atom)
            angle = np.arccos(np.dot(vec12, vec23))

            # dihedral angle between vector(atom1->atom0) and vector(atom2->atom)
            norm0 = np.cross(vec12, vec10)
            norm0 = norm0 / np.sqrt(np.dot(norm0, norm0))
            norm1 = np.cross(vec12, vec23)
            norm1 = norm1 / np.sqrt(np.dot(norm1, norm1))
            #ZULÄSSIG ODER ZU GROB?? WERDEN FALSCHE WINKEL ERMITTELT??
            #MAN KÖNNTE DEN "ÜBERHANG" (ALSO WAS ÜBER 1.0 BZW. -1.0) HINAUSGEHT AUF DEN WINKEL (ALSO 0 BZW. 180 GRAD)
            #IN DIE JEWEILIGE RICHTUNG DRAUFRECHNEN
            temp = np.dot(norm0, norm1)
            if temp > 1.0:
                temp = 1.0
            elif temp < -1.0:
                temp = -1.0
            #dihedral = np.arccos(np.dot(norm0, norm1))
            dihedral = np.arccos(temp)
            if np.dot(np.cross(norm0, norm1), vec12) < 0:
                dihedral = -dihedral

            # add atom with internal coordinates to z-matrix
            new_zmatrix.append([atom3, distance, angle, dihedral])
    

    # BETEILIGTE WBB ATOME MÜSSEN AUS ENTSPRECHENDEN DATENSTRUKTUREN ENTFERNT WERDEN
    def add_monomer_at_acc(self):
        pass
    

    # BETEILIGTE WBB ATOME MÜSSEN AUS ENTSPRECHENDEN DATENSTRUKTUREN ENTFERNT WERDEN
    def add_monomer_at_don(self, cluster: ClusterStructure, atom_to_dock_at: str, docking_atom: str):
        label_shift = len(cluster.coords.keys())
        pos_atom1 = cluster.get_don_vec(atom_to_dock_at)
        pos_atom2 = np.array([])
        pos_atom3 = np.array([])
        hb_vec = pos_atom1 - self.monomer.coords[atom_to_dock_at]
        hb_vec = hb_vec / np.linalg.norm(hb_vec)
        for atom_entry in self.zmatrices[docking_atom]:
            label = atom_entry[0]
            new_label = f"{cluster.get_element(label)}{cluster.get_number(label) + label_shift}"
            new_bond_partners = []
            for bond_partner in self.monomer.bond_partners[label]:
                new_bond_partner_label = f"{cluster.get_element(bond_partner)}{cluster.get_number(bond_partner) + label_shift}"
                new_bond_partners.append(new_bond_partner_label)
            # place first atom of new monomer 
            # using hydrogen bond length and hydrogen bond vector
            if len(atom_entry) == 1:
                new_coords = pos_atom1
            # place second atom of new monomer 
            # using hydrogen bond vector and distance to first atom
            elif len(atom_entry) == 2:
                distance = atom_entry[1]
                pos_atom2 = pos_atom1 + distance * hb_vec
                new_coords = pos_atom2
                vec12 = pos_atom2 - pos_atom1 # PROLBEM, FÜR DIHEDRAL VEC01 BENÖTIGT?? REICHT MAL -1??
                vec12 = vec12 / np.linalg.norm(vec12)
            # place third atom of new monomer 
            # using vector between first and second atom, distance to second atom and angle between 1.-2.-3. atom
            elif len(atom_entry) == 3:
                distance = atom_entry[1]
                angle = atom_entry[2]
                pos_atom3 = pos_atom2 - distance * vec12
                vec23 = pos_atom3 - pos_atom2
                vec23 = vec23 / np.linalg.norm(vec23)
                                                       # DAS ERGIBT NULLVEKTOR!
                #norm_vec = np.cross(-vec12, vec23)    # BELIEBIGEN ANDEREN VEKTOR NEHMEN ODER SMART SO KONSTRUIEREN,
                                                       # DASS MÖGLICHST KEINE CLASHES ENTSTEHEN?
                norm_vec = np.cross(-vec12, pos_atom3) # FÜR DEN MOMENT ERSTMAL IRGENDWAS
                norm_vec = pos_atom2 + norm_vec
                pos_atom3 = rotation(pos_atom3, pos_atom2, norm_vec, angle)
                new_coords = pos_atom3
                vec23 = pos_atom3 - pos_atom2
                vec23 = vec23 / np.linalg.norm(vec23)
            # place fourth, fifth, ..., n-th atom of new monomer
            # using 
            elif len(atom_entry) == 4:
                distance = atom_entry[1]
                angle = atom_entry[2]
                dihedral = atom_entry[3]
                pos_atom4 = pos_atom3 - distance * vec23
                vec34 = pos_atom4 - pos_atom3
                vec34 = vec34 / np.linalg.norm(vec34)
                norm_vec = np.cross(-vec23, pos_atom4)
                norm_vec = pos_atom3 + norm_vec
                pos_atom4 = rotation(pos_atom4, pos_atom3, norm_vec, angle)
                pos_atom4 = rotation(pos_atom4, pos_atom2, pos_atom3, dihedral)
                new_coords = pos_atom4
            cluster.coords[new_label] = new_coords
            cluster.bond_partners[new_label] = new_bond_partners
        self.write_cluster_xyz(cluster.coords)

     
    def adjust_distance(self):
        pass
    def adjust_angle(self):
        pass
    def adjust_dihedral(self):
        pass


    #STATT MONOMER STRUCTURE IN GENERATOR ZU SPEICHERN BEI FUNKTIONEN WIE HIER ÜBERGEBEN
    def old_add_monomer(self, cluster: ClusterStructure, dock_atom: str, docking_atom: str):
        # *** start initialization block ***

        # shift is necessary for shifting number of each label of monomer atoms to new labelnumber
        shift = len(list(cluster.coords))

        #NICHT MEHR NÖTIG MIT update_torsion-FUNKTION
        # after monomer is added the torsion angle will be varied, initialization of new
        # torsion atoms list for this purpose
        #cluster.torsion_atoms = [[]]
        #cluster.torsion_atoms = []

        # position in cluster where new monomer will be added
        position = np.array(cluster.coords[dock_atom])

        # list for coordinates of first atom of z-matrix
        atom0_pos = position
        # list for coordinates of second atom of z-matrix
        atom1_pos = []
        # list for coordinates of third atom of z-matrix
        atom2_pos = []
        # vector first to second atom of z-matrix
        vec10 = []
        # vector second to third atom of z-matrix
        vec12 = []

        # *** end of initialization block ***

        # *** start block for placing every atom of monomer and updating related structural information ***

        # get direction of the hydrogen bond + remove dock atom from list of hydrogen bond atoms (because now occupied)
        # both actions depending on whether dock atom is hydrogen atom donator or acceptor
        #
        #HIER AUC BENDINGS UND BENDING_ATOMS UPDATEN ODER IN SEPARATER FUNKTION??
        if is_hb_acc(dock_atom):
            direction = np.array(cluster.set_acc_vec(dock_atom, cluster.structure[dock_atom][0]))
            cluster.hb_acc.remove(dock_atom)
        else:
            direction = np.array(cluster.set_don_vec(dock_atom))
            cluster.hb_don.remove(dock_atom)

        # loop over every atom of z-matrix of monomer
        for atom in self.zmatrices[docking_atom]:
            # get label of atom and calculate new label with shift
            label = atom[0]
            new_label = get_element(label) + str(shift + get_number(label))

            # add new atom and its neighbors with new labels to structure list of cluster
            cluster.structure[new_label] = []
            #HIER ANPASSEN WENN VERSCHIEDENE MONOMERE VERBAUT WERDEN SOLLEN (self.monomer_structure)
            for neighbor_label in self.monomer_structure[label]:
                new_neighbor_label = get_element(neighbor_label) + str(shift + get_number(neighbor_label))
                cluster.structure[new_label].append(new_neighbor_label)

            # update list of torsion atom and list of hydrogen bond atoms
            # docking atom is not added to list of torsion atoms because it will be part of the rotation axis
            if label != docking_atom:
                #NICHT MEHR NÖTIG MIT update_torsion-FUNKTION
                #cluster.torsion_atoms.append(new_label)
                if is_hb_don(label):
                    cluster.hb_don.append(new_label)
                elif is_hb_acc(label):
                    if is_hb_don(cluster.structure[label][0]):
                        cluster.hb_acc.append(new_label)

            # add first atom (docking atom) of z-matrix to cluster (just based on direction and length of hyrogen bond)
            if len(atom) == 1:
                #ACHTUNG: DAS FUNKTIONIERT SO NUR, WEIL "DIRECTION" BEREITS WBB LÄNGE ENTHÄLT!!
                atom0_pos = position + direction
                cluster.coords[new_label] = list(atom0_pos.copy())

            # add second atom of z-matrix to cluster (based on direction and distance from dock atom)
            elif len(atom) == 2:
                distance = atom[1]
                direction = direction / np.sqrt(np.dot(direction, direction))

                # calculate and add new position of atom
                atom1_pos = atom0_pos + distance * direction
                cluster.coords[new_label] = list(atom1_pos.copy())

                # get vector from first to second atom
                vec10 = atom0_pos - atom1_pos
                vec10 = vec10 / np.sqrt(np.dot(vec10, vec10))

            # add third atom of z-matrix to cluster (based on distance from XYZ and angle XYZ)
            elif len(atom) == 3:
                distance = atom[1]
                angle = atom[2]

                # calculate and add new position of atom
                atom2_pos = atom1_pos + distance * vec10
                norm = np.cross(vec10, atom2_pos)
                norm = atom1_pos + norm
                atom2_pos = rotation(norm, atom1_pos, angle, atom2_pos, True)
                cluster.coords[new_label] = list(atom2_pos.copy())

                # get vector from second to third atom
                vec12 = atom2_pos - atom1_pos
                vec12 = vec12 / np.sqrt(np.dot(vec12, vec12))

            # add every further atom after first three atoms of z-matrix (based on distance from XYZ,
            # angle XYZ and dihedral angle XYZ)
            elif len(atom) == 4:
                distance = atom[1]
                angle = atom[2]
                dihedral = atom[3]

                # calculate and add position of atom
                new_atom_pos = atom2_pos + distance * vec12
                norm0 = np.cross(new_atom_pos, vec10)
                norm0 = atom2_pos + norm0
                new_atom_pos = rotation(norm0, atom2_pos, angle, new_atom_pos, True)
                new_atom_pos = rotation(atom2_pos, atom1_pos, dihedral, new_atom_pos, True)
                cluster.coords[new_label] = list(new_atom_pos)

        # get new lagel of docking atom
        docking_atom_new_label = get_element(docking_atom) + str(get_number(docking_atom) + shift)

        #NICHT MEHR NÖTIG MIT update_torsion-FUNKTION
        # update rotation axis for torsion angle variation
        #cluster.torsions = [dock_atom, docking_atom_new_label]

        # add docking atom to neighbor list of dock atom
        cluster.structure[dock_atom].append(docking_atom_new_label)

        # add dock atom to neighbor list of docking atom
        cluster.structure[docking_atom_new_label].append(dock_atom)

        # update number of monomers building up the current cluster
        cluster.monomers += 1


    def update_torsion(self, cluster: ClusterStructure, dock_atom: str, docking_atom: str):
        cluster.torsions = [dock_atom, docking_atom]
        cluster.torsion_atoms = [atom[0] for atom in self.zmatrices[docking_atom]]


    #def update_bending(self, cluster: ClusterStructure, dock_atom: str, docking_atom: str):
    def update_bending(self, cluster: ClusterStructure, don: str, acc: str):
        neighbor = cluster.structure[don][0]
        self.set_bendings(cluster, don, acc, neighbor)
        status = {atom: "UNKNOWN" for atom in cluster.structure}
        self.set_bending_atoms(cluster, don, status, acc)
        #if is_hb_don(dock_atom):
        #    neighbor = cluster.structure[dock_atom][0]
        #    self.set_bendings(cluster, dock_atom, docking_atom, neighbor)
        #    status = {atom: "UNKNOWN" for atom in cluster.structure}
        #    self.set_bending_atoms(cluster, dock_atom, status, docking_atom)
        #else:
        #    neighbor = cluster.structure[docking_atom][0]
        #    self.set_bendings(cluster, docking_atom, dock_atom, neighbor)
        #    status = {atom: "UNKNOWN" for atom in cluster.structure}
        #    self.set_bending_atoms(cluster, docking_atom, status, dock_atom)


    #WIE SOLL MIT DONATOREN MIT EINEM EINZIGEN NACHBARN VERFAHREN WERDEN (Z.B. CARBOXYLSAUERSTOFF)?? DAFÜR
    #SEPZIALFALL ANLEGEN??
    def set_bendings(self, cluster: ClusterStructure, don: str, acc: str=None, neighbor: str = None):
        index = 0
        # get label of neighbor atom of donator atom if none is given
        if not neighbor:
            # every neighbor could be used, take first one in structure list
            neighbor = cluster.structure[don][index]
        # neighbor should be different from acceptor atom, if this is not the case for first neighbor
        # in structure list take the next (second) one
        while neighbor == acc:
            index += 1
            neighbor = cluster.structure[don][index]

        # number of neighbor atoms of donator atom
        n = len(cluster.structure[don])
        # store coordinates of donator atom in don variable
        don = np.array(cluster.coords[don])
        # store coordinates of acceptor atom in acc variable
        acc = np.array(cluster.coords[acc])
        # store coordinats of neighbor atom in neighbor variable
        neighbor = np.array(cluster.coords[neighbor])

        # get vector from donator atom to acceptor atom
        vec_a = acc - don
        vec_a = vec_a / np.sqrt(np.dot(vec_a, vec_a))
        # get vector from donator atom to neighbor atom
        vec_b = neighbor - don
        vec_b = vec_b / np.sqrt(np.dot(vec_b, vec_b))

        # get vector which is orthogonal to vec_a and vec_b and normalize it
        axis = np.cross(vec_a, vec_b)
        axis = axis / np.sqrt(np.dot(axis, axis))
        # get first rotation axis by rotating orthogonal vector around vec_a
        #axis = rotation(acc, don, np.deg2rad(-1 * (n+1)/n * 180), don + axis)
        axis = rotation(acc, don, 180/n, don + axis)
        # calculate rotation axes (position of axis atoms) by rotating axis around vec_a and store them in bendings
        for i in range(n):
            cluster.bendings.append([list(rotation(acc, don, 360/n * i, axis)), list(don)])


    def set_bending_atoms(self, cluster: ClusterStructure, current_atom: str, status: dict, acc: str):
        status[current_atom] = "KNOWN"
        for atom in cluster.structure[current_atom]:
            if status[atom] == "UNKNOWN":
                if atom != acc:
                    cluster.bending_atoms.append(atom)
                    self.set_bending_atoms(cluster, atom, status, acc)


    def conformers(self, oligomer: ClusterStructure, monomer: ClusterStructure, max_size: int):
        # base case, Torsions- und Biegewinkel für jede WBB wurde berechnet
        if oligomer.monomers == max_size:
            self.output(oligomer.coords)
            self.counter += 1
            return
        # es wurden noch nicht alle Torsionswinkel berechnet, Bindung index in torsions ist an der Reihe
        else:
            """
            #NOCHMAL DIESE SCHLEIFENKOMBI, DA KOMBIS AUCH ANDERSRUM DURCHLAUFEN WERDEN MÜSSEN!!!
            #DANN AUS DARAUF FOLGENDEM CODEBLOCK EIGENE FUNKTION MACHEN, DA SONST QUASI DERSELBE CODE
            #ZWEIMAL VORKÄME (BLOCK NACH UPDATE FUNKTIONEN)
            #DANN ALLERDINGS FÜR DIESEN TEIL DOPPELT (?) SO VIELE FUNKTIONSAUFRUFE WIE ZUVOR, WIE GROß IST
            #GEFAHR EINES STACK OVERFLOWS??
            for don in oligomer.hb_don:
                for acc in monomer.hb_acc:
                    NewCluster = deepcopy(oligomer)
                    self.add_monomer(NewCluster, don, acc)
                    self.update_torsion(NewCluster, don, acc)
                    self.update_bending(NewCluster, don, acc)

                    # vectors forming rotation axis for torsion angle
                    #vec11 = oligomer.coords[oligomer.torsions[0]].copy()
                    #vec12 = oligomer.coords[oligomer.torsions[1]].copy()
                    vec11 = NewCluster.coords[NewCluster.torsions[0]].copy()
                    vec12 = NewCluster.coords[NewCluster.torsions[1]].copy()

                    #WINKEL UNGESCHICKT GEWÄHLT WENN ZB HYDRONIUMIONEN BETRACHTET WERDEN
                    #WAS WENN ROTATION GAR NICHT SINNVOLL IST ZB HYDROXIDIONEN
                    for torsion_angle in [0, 120, 240]:
                        for atom in NewCluster.torsion_atoms:
                            coords = NewCluster.coords[atom].copy()
                            coords = rotation(vec12, vec11, torsion_angle, coords)
                            NewCluster.coords[atom] = list(coords.copy())
                        for bending in NewCluster.bendings:
                            vec21 = bending[0].copy()
                            vec22 = bending[1].copy()
                            for bending_angle in [0, 45]:
                                for atom in NewCluster.bending_atoms:
                                    coords = NewCluster.coords[atom].copy()
                                    coords = rotation(vec22, vec21, bending_angle, coords)
                                    NewCluster.coords[atom] = list(coords.copy())

                                #if not NewCluster.clashes():
                                #    self.conformers(NewCluster, monomer, max_size)
                                self.conformers(NewCluster, monomer, max_size)
            """
            for acc in oligomer.hb_acc:
                for don in monomer.hb_don:
                    NewCluster = deepcopy(oligomer)
                    self.add_monomer(NewCluster, acc, don)
                    self.update_torsion(NewCluster, acc, don)
                    self.update_bending(NewCluster, don, acc)
                    #temp = NewCluster.coords[acc]
                    #NewCluster.coords["Br" + str(get_number(acc))] = temp
                    #temp = NewCluster.coords[don]
                    #NewCluster.coords["Br" + str(get_number(don))] = temp

                    # vectors forming rotation axis for torsion angle
                    # vec11 = oligomer.coords[oligomer.torsions[0]].copy()
                    # vec12 = oligomer.coords[oligomer.torsions[1]].copy()
                    vec11 = NewCluster.coords[NewCluster.torsions[0]].copy()
                    vec12 = NewCluster.coords[NewCluster.torsions[1]].copy()

                    # WINKEL UNGESCHICKT GEWÄHLT WENN ZB HYDRONIUMIONEN BETRACHTET WERDEN
                    # WAS WENN ROTATION GAR NICHT SINNVOLL IST ZB HYDROXIDIONEN
                    for torsion_angle in [0, 120, 240]:
                        for atom in NewCluster.torsion_atoms:
                            coords = NewCluster.coords[atom].copy()
                            coords = rotation(vec12, vec11, torsion_angle, coords)
                            NewCluster.coords[atom] = list(coords.copy())
                        for bending in NewCluster.bendings:
                            vec21 = bending[0].copy()
                            vec22 = bending[1].copy()
                            for bending_angle in [45]:
                                for atom in NewCluster.bending_atoms:
                                    coords = NewCluster.coords[atom].copy()
                                    coords = rotation(vec22, vec21, bending_angle, coords)
                                    #NewCluster.coords[atom] = list(coords.copy())
                                    NewCluster.coords["Cl" + str(get_number(atom))] = list(coords.copy())

                                #if not NewCluster.clashes():
                                #    self.conformers(NewCluster, monomer, max_size)
                                self.conformers(NewCluster, monomer, max_size)


    def write_cluster_xyz(self, cluster_coords: dict) -> None:
        with open("out.xyz", "w") as xyzfile:
            n_atoms = len(cluster_coords.keys())
            print(n_atoms, file=xyzfile, end="\n\n")
            for atom, coords in cluster_coords.items():
                element = self.monomer.get_element(atom)
                x = coords[0]
                y = coords[1]
                z = coords[2]
                print(f"{element}   {x:.5f}   {y:.5f}   {z:.5f}", file=xyzfile)


    def output(self, coords, counter: int = None):
        if not counter:
            counter = self.counter
        if not os.path.exists("../Output"):
            os.makedirs("../Output")
        with open("../Output/conformer" + str(counter) + ".xyz", "w") as outfile:
            outfile.write(str(len(coords)))
            outfile.write("\n\n")
            for atom in coords:
                outfile.write(get_element(atom))
                outfile.write("\t")
                print("{:18.15f} \t {:18.15f} \t {:18.15f}".format(coords[atom][0], coords[atom][1],
                                                                   coords[atom][2]), file=outfile)
