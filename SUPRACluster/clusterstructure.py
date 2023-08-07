from SUPRAConformer.structure import Structure
from SUPRACluster.hydrogen_bond_donator import HydrogenBondDonator
from SUPRACluster.hydrogen_bond_acceptor import HydrogenBondAcceptor
from utils.helper import get_element



class ClusterStructure(Structure):

    def __init__(self, file: str=None):
        super(ClusterStructure, self).__init__(file)

        # number of monomers that the cluster consists of
        self.monomers = 1

        # hydrogen bond donators which are free to bind to a hydrogen bond acceptor
        #SHOW GENERAL STRUCTURES
        self.hydrogen_bond_acceptors = {}

        # hydrogen bond acceptors which are free to bind to a hydrogen bond donator
        #SHOW GENERAL STRUCTURES
        self.hydrogen_bond_donators = {}

        self.hydrogen_bonds = []

        # NOCH (WO)ANDERS UNTERBRINGEN UND AUF ANDERE ELEMENTE ANPASSEN
        self.hb_len = 1.1

        self.find_hb_accs_dons()

    

    # find and store atoms that form hydrogen bonds (hbs), so hb donors or acceptors
    def find_hb_accs_dons(self) -> None:
        for atom in self.coords.keys():
            if get_element(atom) in ["O", "N"]:
                new_don = HydrogenBondDonator(atom, self.coords[atom])
                new_don.get_don_vec(self.bond_partners[atom], self.coords)
                self.hydrogen_bond_donators[atom] = new_don
                for neighbor in self.bond_partners[atom]:
                    if get_element(neighbor) == "H":
                        new_acc = HydrogenBondAcceptor(neighbor, self.coords[neighbor])
                        new_acc.get_acc_vec(self.coords[self.bond_partners[neighbor][0]])
                        self.hydrogen_bond_acceptors[neighbor] = new_acc
