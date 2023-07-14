from utils.bond import Bond



class Torsion(Bond):

    def __init__(self) -> None:
        super(Bond, self).__init__()
        self.torsion_atoms = []
        self.rot_atoms1 = []
        self.rot_atoms2 = []
        self.sym_rot_atoms1 = []
        self.sym_rot_atoms2 = []
        self.rot_sym1 = 1
        self.rot_sym2 = 1
        self.rot_angles = []