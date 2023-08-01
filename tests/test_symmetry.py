import os

import numpy as np

from utils.helper import get_number, geometric_center
from utils.symmetry import Symmetry
from utils.rotationaxis import RotationAxis
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
    sym = Symmetry()
    mol1 = Structure(os.path.join(files, "Alanin.xyz"))
    mol2 = Structure(os.path.join(files, "Benzol.xyz"))

    atoms = ["C2", "H3", "H4", "H5"]
    from_coords = mol1.coords["C0"]
    to_coords = mol1.coords["C2"]
    rotaxis1 = RotationAxis(from_coords, to_coords)

    check1 = sym.rot_sym_along_bond(mol1, rotaxis1, atoms, 3)
    check2 = sym.rot_sym_along_bond(mol1, rotaxis1, atoms, 2)

    atoms = ["N6", "H8", "H9"]
    from_coords = mol1.coords["C0"]
    to_coords = mol1.coords["N6"]
    rotaxis2 = RotationAxis(from_coords, to_coords)

    check3 = sym.rot_sym_along_bond(mol1, rotaxis2, atoms, 1)
    check4 = sym.rot_sym_along_bond(mol1, rotaxis2, atoms, 2)

    atoms = list(mol2.coords.keys())
    from_coords = geometric_center(mol2.coords)
    vec1 = (from_coords - mol2.coords["C0"])
    vec1 = vec1 / np.linalg.norm(vec1)
    vec2 = (from_coords - mol2.coords["C1"])
    vec2 = vec2 / np.linalg.norm(vec2)
    norm_vec = np.cross(vec1, vec2)
    norm_vec = norm_vec / np.linalg.norm(norm_vec)
    to_coords = from_coords + norm_vec
    rotaxis3 = RotationAxis(from_coords, to_coords)

    check5 = sym.rot_sym_along_bond(mol2, rotaxis3, atoms, 6)
    check6 = sym.rot_sym_along_bond(mol2, rotaxis3, atoms, 3)

    assert (check1 and not check2 and check3 and not check4 and check5 and check6)



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
