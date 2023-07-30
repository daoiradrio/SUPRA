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
    gen = ConformerGenerator()
    mol1 = Structure(os.path.join(files, "Alanin.xyz"))
    mol2 = Structure(os.path.join(files, "Tyrosin.xyz"))

    gen._find_torsions(mol1.bonds, mol1.bond_partners)

    central_bond = gen.central_torsions[0]
    methyl_bond = gen.methyl_torsions[0]

    check1 = False
    for bond in gen.central_torsions:
        check1 = (bond.atom1 == "C0" and bond.atom2 == "C1")
        if check1:
            break
    check2 = False
    for bond in gen.methyl_torsions:
        check2 = (bond.atom1 == "C0" and bond.atom2 == "C2")
        if check2:
            break

    gen._find_torsions(mol2.bonds, mol2.bond_partners)

    check3 = False
    for bond in gen.terminal_torsions:
        check3 = (bond.atom1 == "C1" and bond.atom2 == "O6")
        if check3:
            break

    assert (check1 and check2 and check3)



def test_find_cylces():
    gen = ConformerGenerator()
    mol = Structure(os.path.join(files, "Propylcyclohexan.xyz"))

    gen._find_torsions(mol.bonds, mol.bond_partners)

    check1 = False
    for bond in gen.central_torsions:
        check1 = (bond.atom1 == "C8" and bond.atom2 == "C16")
        if check1:
            break
    check2 = False
    for bond in gen.central_torsions:
        check2 = (bond.atom1 == "C16" and bond.atom2 == "C18")
        if check2:
            break
    check3 = False
    for bond in gen.central_torsions:
        check3 = (bond.atom1 == "C0" and bond.atom2 == "C1")
        if check3:
            break

    gen._find_cycles(mol.bond_partners)

    check4 = False
    for bond in gen.central_torsions:
        check4 = (bond.atom1 == "C8" and bond.atom2 == "C16")
        if check4:
            break
    check5 = False
    for bond in gen.central_torsions:
        check5 = (bond.atom1 == "C16" and bond.atom2 == "C18")
        if check5:
            break
    check6 = False
    for bond in gen.central_torsions:
        check6 = (bond.atom1 == "C0" and bond.atom2 == "C1")
        if check6:
            break

    assert (check1 and check2 and check3 and check4 and check5 and not check6)



def test_find_peptidebonds():
    gen = ConformerGenerator()
    mol = Structure(os.path.join(files, "Tyr-Ala-Trp.xyz"))

    gen._find_torsions(mol.bonds, mol.bond_partners)

    check1 = False
    for bond in gen.central_torsions:
        check1 = (bond.atom1 == "C14" and bond.atom2 == "N22")
        if check1:
            break
    check2 = False
    for bond in gen.central_torsions:
        check2 = (bond.atom1 == "C29" and bond.atom2 == "N32")
        if check2:
            break
    
    gen._find_peptidebonds(mol.coords, mol.bond_partners)

    check3 = False
    for bond in gen.central_torsions:
        check3 = (bond.atom1 == "C14" and bond.atom2 == "N22")
        if check3:
            break
    check4 = False
    for bond in gen.central_torsions:
        check4 = (bond.atom1 == "C29" and bond.atom2 == "N32")
        if check4:
            break
    
    assert (check1 and check2 and not check3 and not check4)



def test_find_clashes():
    gen = ConformerGenerator()
    mol1 = Structure(os.path.join(files, "Tyrosin.xyz"))
    mol2 = Structure(os.path.join(files, "Tyrosin_clash.xyz"))

    check1 = gen._find_clashes(mol1.bond_partners, mol1.coords)
    check2 = gen._find_clashes(mol1.bond_partners, mol2.coords)

    assert (not check1 and check2)
