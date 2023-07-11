import numpy as np

from SUPRAConformer.structure import Structure
from utils.hydrogen_bond_donator import HydrogenBondDonator
from utils.hydrogen_bond_acceptor import HydrogenBondAcceptor
from utils.helper import get_element



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
        self.hydrogen_bond_acceptors = []

        # hydrogen bond acceptors which are free to bind to a hydrogen bond donator
        #SHOW GENERAL STRUCTURES
        self.hydrogen_bond_donators = []

        # NOCH (WO)ANDERS UNTERBRINGEN UND AUF ANDERE ELEMENTE ANPASSEN
        self.hb_len = 1.1

        self.find_hb_accs_dons()

    

    # find and store atoms that form hydrogen bonds (hbs), so hb donors or acceptors
    def find_hb_accs_dons(self) -> None:
        for atom in self.coords.keys():
            if get_element(atom) in ["O", "N"]:
                new_don = HydrogenBondDonator(atom, self.coords[atom])
                new_don.get_don_vec(self.bond_partners[atom], self.coords)
                self.hydrogen_bond_donators.append(new_don)
                for neighbor in self.bond_partners[atom]:
                    if get_element(neighbor) == "H":
                        new_acc = HydrogenBondAcceptor(neighbor, self.coords[neighbor])
                        new_acc.get_acc_vec(self.coords[self.bond_partners[neighbor][0]])
                        self.hydrogen_bond_acceptors.append(new_acc)
