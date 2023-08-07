class Bond:

    def __init__(self, atom1: str=None, atom2: str=None, bond_order: int=0) -> None:
        self.bond_order = bond_order
        self.atom1 = atom1
        self.atom2 = atom2
