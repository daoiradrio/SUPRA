import os

from SUPRAConformer.structure import Structure
from utils.helper import get_element



files = "testcases/"
try:
    f = open(os.path.join(files, "Alanin.xyz"))
    f.close()
except:
    files = "tests/testcases/"



def test_read_xyz():
    mol = Structure(os.path.join(files, "Alanin.xyz"))
    atoms = [get_element(atom) for atom in mol.coords.keys()]
    C_count = atoms.count("C")
    H_count = atoms.count("H")
    O_count = atoms.count("O")
    N_count = atoms.count("N")
    assert (
        mol.number_of_atoms == 13 and 
        C_count == 3 and
        H_count == 7 and
        O_count == 2 and
        N_count == 1
    )



def test_get_connectivity():
    mol = Structure(os.path.join(files, "Tyrosin.xyz"))
    atoms = {
        atom: [
            get_element(bond_partner) for bond_partner in mol.bond_partners[atom]
        ] for atom in mol.bond_partners.keys()
    }
    C0 = (
        atoms["C0"].count("C") == 2 and 
        atoms["C0"].count("H") == 1
    )
    C1 = (
        atoms["C1"].count("C") == 2 and 
        atoms["C1"].count("O") == 1
    )
    C2 = (
        atoms["C2"].count("C") == 2 and
        atoms["C2"].count("H") == 1 
    )
    C3 = (
        atoms["C3"].count("C") == 2 and
        atoms["C3"].count("H") == 1
    )
    C4 = (
        atoms["C4"].count("C") == 3
    )
    C5 = (
        atoms["C5"].count("C") == 2 and
        atoms["C5"].count("H") == 1
    )
    O6 = (
        atoms["O6"].count("C") == 1 and
        atoms["O6"].count("H") == 1
    )
    C12 = (
        atoms["C12"].count("C") == 2 and
        atoms["C12"].count("H") == 2
    )
    C13 = (
        atoms["C13"].count("C") == 2 and
        atoms["C13"].count("H") == 1 and
        atoms["C13"].count("N") == 1
    )
    C16 = (
        atoms["C16"].count("C") == 1 and
        atoms["C16"].count("O") == 2
    )
    O17 = (
        atoms["O17"].count("C") == 1
    )
    O18 = (
        atoms["O18"].count("C") == 1 and
        atoms["O18"].count("H") == 1
    )
    N20 = (
        atoms["N20"].count("C") == 1 and
        atoms["N20"].count("H") == 2
    )
    assert (
        C0 and C1 and C2 and C3 and C4 and C5 and O6 and 
        C12 and C13 and C16 and O17 and O18 and N20
    )