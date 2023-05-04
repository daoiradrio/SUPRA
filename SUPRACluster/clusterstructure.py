import numpy as np

from SUPRAConformer.structure import Structure
from utils.helper import is_hb_don, is_hb_acc



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

        # number of monomers that the cluster consists of
        self.monomers = 1

        # hydrogen bond donators which are free to bind to a hydrogen bond acceptor
        #SHOW GENERAL STRUCTURES
        self.hb_don = list()

        # hydrogen bond acceptors which are free to bind to a hydrogen bond donator
        #SHOW GENERAL STRUCTURES
        self.hb_acc = list()

        # NOCH (WO)ANDERS UNTERBRINGEN UND AUF ANDERE ELEMENTE ANPASSEN
        self.hb_len = 1.1

        self.find_hbs()

        #DATENSTRUKTUR UM WBB (ATOMLABEL UND POSITION) ZU SPEICHERN UM ANHAND DESSEN KONFORMERE ERZEUGEN ZU KÖNNEN

    
    # find and store atoms that form hydrogen bonds (hbs), so hb donors or acceptors
    def find_hbs(self) -> None:
        for atom in self.coords.keys():
            if self.get_element(atom) in ["O", "N"]:
                self.hb_don.append(atom)
                #self.hb_don_vec[atom] = self.get_don_vec(atom)
                for neighbor in self.bond_partners[atom]:
                    if self.get_element(neighbor) == "H":
                        self.hb_acc.append(neighbor)
                        #self.hb_acc_vec[neighbor] = self.get_acc_vec(neighbor)

    
    def get_don_vec(self, don: str) -> np.array:
        vecs = []
        for neighbor in self.bond_partners[don]:
            new_vec = self.coords[don] - self.coords[neighbor]
            new_vec = new_vec / np.linalg.norm(new_vec)
            vecs.append(new_vec)
        new_hb = self.coords[don].copy()
        for vec in vecs:
            new_hb += vec
        return new_hb


    def get_acc_vec(self, acc: str) -> np.array:
        neighbor = self.bond_partners[acc][0]
        new_hb = self.coords[acc] - self.coords[neighbor]
        new_hb = new_hb / np.linalg.norm(new_hb)
        new_hb = self.coords[acc] + new_hb
        return new_hb
