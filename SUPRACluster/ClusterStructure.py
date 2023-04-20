from Structure import Structure
import numpy as np
from Helper import get_number, is_hb_don, is_hb_acc



class ClusterStructure(Structure):

    def __init__(self, file: str=None):
        super(ClusterStructure, self).__init__(file)

        """
        *** INHERITED FROM STRUCTURE-CLASS ***
        
        coordinates of all cluster atoms
        coords = {E0:[x0, y0, z0], E1:[x1, y1, z1], ... }
        self.coords = dict()
        
        connectivity of all atoms in input structure (undirected graph)
        structure = {E0:[E1, E2, E3, ... ], E1:[E0, E4, ], ... }
        self.structure = dict()
        
        # nonterminal single bonds, unfiltered rotatable bonds
        # torsions = [[E0, E1], [E2, E3], ... ]
        self.torsions = list()
        
        self.torsion_atoms = list()
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


    def find_hbs(self):
        for atom in self.coords:
            if is_hb_don(atom):
                self.hb_don.append(atom)
                self.hb_don_vec[get_number(atom)] = self.set_don_vec(atom)
                for neighbor in self.structure[atom]:
                    if is_hb_acc(neighbor):
                        self.hb_acc.append(neighbor)
                        self.hb_acc_vec[get_number(neighbor)] = self.set_acc_vec(neighbor, self.structure[neighbor][0])


    def set_don_vec(self, don: str) -> list:
        vecs = list()
        don_coord = np.array(self.coords[don])
        for neighbor in self.structure[don]:
            neighbor_coord = np.array(self.coords[neighbor])

            new_vec = don_coord - neighbor_coord
            len = np.sqrt(np.dot(new_vec, new_vec))
            new_vec = new_vec / len

            vecs.append(new_vec)

        new_coord = np.array(self.coords[don])

        for vec in vecs:
            new_coord += vec

        new_hb = new_coord - don_coord
        new_hb_len = np.sqrt(np.dot(new_hb, new_hb))

        #HIER ELEMENTSPEZFISCHE LÄNGE VON WBB VERWENDEN --> QUELLE SUCHEN
        scaling_factor = self.hb_len / new_hb_len

        new_hb = new_hb * scaling_factor

        return list(new_hb)


    def set_acc_vec(self, acc: str, neighbor: str) -> list:
        vec2 = np.array(self.coords[acc])
        vec1 = np.array(self.coords[neighbor])
        direction = vec2 - vec1
        len = np.sqrt(np.dot(direction, direction))
        scaling_factor = self.hb_len / len
        direction = direction * scaling_factor

        return list(direction)
