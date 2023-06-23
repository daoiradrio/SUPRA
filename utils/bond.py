class Bond:

    def __init__(self) -> None:
        self.bond_order = 0
        self.atom1 = None
        self.atom2 = None
        self.torsion_atoms = []
        self.rot_atoms1 = []
        self.rot_atoms2 = []
        self.sym_rot_atoms1 = []
        self.sym_rot_atoms2 = []
        self.rot_sym1 = 1
        self.rot_sym2 = 1
        self.rot_angles = []