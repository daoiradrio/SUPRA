import os

from SUPRAConformer.structure import Structure
from SUPRAConformer.conformergenerator import ConformerGenerator



files = "testcases/"
try:
    f = open(os.path.join(files, "Alanin.xyz"))
    f.close()
except:
    files = "tests/testcases/"



def test_find_torsions():
    mol1 = Structure(os.path.join(files, "Alanin.xyz"))
    mol2 = Structure(os.path.join(files, "Tyrosin.xyz"))
    gen = ConformerGenerator()

    gen._find_torsions(mol1.bonds, mol1.bond_partners)

    central_bond = gen.central_torsions[0]
    methyl_bond = gen.methyl_torsions[0]

    gen._find_torsions(mol2.bonds, mol2.bond_partners)

    terminal_bond1 = gen.terminal_torsions[0]
    terminal_bond2 = gen.terminal_torsions[1]
    terminal_bond3 = gen.terminal_torsions[2]

    check1 = (central_bond.atom1 == "C0" and central_bond.atom2 == "C1")
    check2 = (methyl_bond.atom1 == "C0" and methyl_bond.atom2 == "C2")
    check3 = (terminal_bond1.atom1 ==  "C1" and terminal_bond1.atom2 ==  "O6")
    check4 = (terminal_bond2.atom1 == "C13" and terminal_bond2.atom2 == "N20")
    check5 = (terminal_bond3.atom1 == "C16" and terminal_bond3.atom2 == "O18")

    assert (check1 and check2 and check3 and check4 and check5)



def test_find_cylces():
    assert True

def test_peptidebonds():
    assert True

def test_clashes():
    assert True