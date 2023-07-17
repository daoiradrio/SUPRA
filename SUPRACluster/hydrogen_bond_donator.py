import numpy as np



class HydrogenBondDonator:

    def __init__(self, atom_label: str, atom_coords: np.array) -> None:
        self.atom_label = atom_label
        self.atom_coords = atom_coords
        self.vec = None
    

    
    def get_don_vec(self, bond_partners: list, coords: dict):
        vecs = []
        for neighbor in bond_partners:
            new_vec = self.atom_coords - coords[neighbor]
            new_vec = new_vec / np.linalg.norm(new_vec)
            vecs.append(new_vec)
        self.vec = np.array([0.0, 0.0, 0.0])
        for v in vecs:
            self.vec += v
