import numpy as np



class HydrogenBondAcceptor:

    def __init__(self, atom_label: str, atom_coords: np.array) -> None:
        self.atom_label = atom_label
        self.atom_coords = atom_coords
        self.vec = None
    


    def get_acc_vec(self, bond_partner_coords: np.array):
        new_hb = self.atom_coords - bond_partner_coords
        new_hb = new_hb / np.linalg.norm(new_hb)
        self.vec = self.atom_coords + new_hb
