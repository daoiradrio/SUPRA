import numpy as np

from SUPRAConformer.structure import Structure
from utils.helper import get_number, is_hb_don, is_hb_acc



class ClusterStructure(Structure):

    def __init__(self, file: str=None):
        super(ClusterStructure, self).__init__(file)

        """
        *** INHERITED FROM STRUCTURE-CLASS ***
        
        numer of atoms
        self.number_of_atoms: int

        coordinates of all atoms
        self.coords = {A0:[x0, y0, z0], A1:[x1, y1, z1], ... }
        
        connectivity of all atoms (undirected graph)
        self.bond_partners = {A0:[A1, A2, A3, ... ], A1:[A0, A4, ], ... }
        
        all covalent bonds (no doubles with just switched atom order)
        self.bonds = [[A0, A1], [A0, A2], [A1, A4], ... ]
        
        bond orders of all covalent bonds
        self.bond_orders = {[A0, A1]: 1, [A1, A4]: 2, ... }
        """

        self.bendings = list()
        self.bending_atoms = list()

        # number of monomers that the cluster consists of
        self.monomers = 1

        # hydrogen bond donators which are free to bind to a hydrogen bond acceptor
        #SHOW GENERAL STRUCTURES
        self.hb_don = list()
        self.hb_don_vec = dict()

        # hydrogen bond acceptors which are free to bind to a hydrogen bond donator
        #SHOW GENERAL STRUCTURES
        self.hb_acc = list()
        self.hb_acc_vec = dict()

        # NOCH (WO)ANDERS UNTERBRINGEN UND AUF ANDERE ELEMENTE ANPASSEN
        self.hb_len = 1.1

        self.find_hbs()

        #DATENSTRUKTUR UM WBB (ATOMLABEL UND POSITION) ZU SPEICHERN UM ANHAND DESSEN KONFORMERE ERZEUGEN ZU KÖNNEN

    
    # find and store atoms that form hydrogen bonds (hbs), so hb donors or acceptors
    def find_hbs(self):
        #for atom in self.coords.keys():
        #    if get_element(atom) == "H":
        #        pass

        for atom in self.coords:
            if is_hb_don(atom):
                self.hb_don.append(atom)
                self.hb_don_vec[get_number(atom)] = self.set_don_vec(atom)
                for neighbor in self.bond_partners[atom]:
                    if is_hb_acc(neighbor):
                        self.hb_acc.append(neighbor)
                        self.hb_acc_vec[get_number(neighbor)] = self.set_acc_vec(neighbor, self.bond_partners[neighbor][0])


    def set_don_vec(self, don: str) -> list:
        vecs = []
        don_coord = self.coords[don]

        for neighbor in self.bond_partners[don]:
            neighbor_coord = self.coords[neighbor]
            new_vec = don_coord - neighbor_coord
            len_new_vec = np.sqrt(np.dot(new_vec, new_vec))
            new_vec = new_vec / len_new_vec
            vecs.append(new_vec)

        new_coord = self.coords[don].copy()
        for vec in vecs:
            new_coord += vec

        new_hb = new_coord - don_coord
        new_hb_len = np.sqrt(np.dot(new_hb, new_hb))
        #HIER ELEMENTSPEZFISCHE LÄNGE VON WBB VERWENDEN --> QUELLE SUCHEN
        scaling_factor = self.hb_len / new_hb_len
        new_hb = new_hb * scaling_factor

        return new_hb


    def set_acc_vec(self, acc: str, neighbor: str) -> list:
        vec2 = self.coords[acc]
        vec1 = self.coords[neighbor]
        new_hb = vec2 - vec1
        len_bond = np.sqrt(np.dot(new_hb, new_hb))
        scaling_factor = self.hb_len / len_bond
        new_hb = new_hb * scaling_factor
        return new_hb
