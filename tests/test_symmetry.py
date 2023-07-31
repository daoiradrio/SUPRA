import os

from utils.helper import get_number
from utils.symmetry import Symmetry
from SUPRAConformer.structure import Structure
from SUPRAConformer.conformergenerator import ConformerGenerator



files = "testcases/"
try:
    f = open(os.path.join(files, "Alanin.xyz"))
    f.close()
except:
    files = "tests/testcases/"



def test_find_torsion_group():
    sym = Symmetry()
    gen = ConformerGenerator()
    mol1 = Structure(os.path.join(files, "para-Phenylphenol.xyz"))
    mol2 = Structure(os.path.join(files, "Tyrosin.xyz"))

    gen._find_torsions(mol1.bonds, mol1.bond_partners)

    status = {atom: "UNKNOWN" for atom in mol1.coords.keys()}
    atoms = []
    sym._find_torsion_group(mol1.bond_partners, gen.torsions, "C2", "O10", status, atoms)

    atoms = sorted(atoms, key=lambda x: get_number(x))
    check1 = (atoms == [
        'C0', 'C1', 'C2', 'C3', 'C4', 
        'C5', 'C5', 'H6', 'H7', 'H8', 
        'H9', 'C12', 'C12', 'C13', 'C14', 
        'C15', 'C16', 'C17', 'C17', 'H18', 
        'H19', 'H20', 'H21', 'H22'
        ])

    gen._find_torsions(mol2.bonds, mol2.bond_partners)

    status = {atom: "UNKNOWN" for atom in mol2.coords.keys()}
    atoms = []
    sym._find_torsion_group(mol2.bond_partners, gen.torsions, "C16", "C13", status, atoms)

    atoms = sorted(atoms, key=lambda x: get_number(x))
    check2 = (atoms == ["O17", "O18", "H19"])

    assert (check1 and check2)



def test_rot_sym_along_bond():
    assert True



def test_rot_order_along_bond():
    sym = Symmetry()
    mol1 = Structure(os.path.join(files, "Alanin.xyz"))
    mol2 = Structure(os.path.join(files, "Benzol.xyz"))

    atoms = ["C2", "H3", "H4", "H5"]
    from_coords = mol1.coords["C0"]
    to_coords = mol1.coords["C2"]

    rot_order1 = sym.rot_order_along_bond(mol1, atoms, from_coords, to_coords)

    atoms = list(mol2.coords.keys())
    from_coords = mol2.coords["H8"]
    to_coords = mol2.coords["C0"]

    rot_order2 = sym.rot_order_along_bond(mol2, atoms, from_coords, to_coords)

    atoms = ["N6", "H8", "H9"]
    from_coords = mol1.coords["C0"]
    to_coords = mol1.coords["N6"]

    rot_order3 = sym.rot_order_along_bond(mol1, atoms, from_coords, to_coords)
    
    assert (rot_order1 == 3 and rot_order2 == 2 and rot_order3 == 1)
