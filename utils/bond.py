class Bond:

    def __init__(self, atom1: str, atom2: str, bond_order: int) -> None:
        self.bond_order = bond_order
        self.atom1 = atom1
        self.atom2 = atom2
