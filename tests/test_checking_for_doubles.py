import os
from SUPRAConformer.structure import Structure
from utils.analyzer import Analyzer



files = "testcases/"
try:
    f = open(os.path.join(files, "Alanin.xyz"))
    f.close()
except:
    files = "tests/testcases/"



def test_loose_doubles_check():
    mol1 = Structure(os.path.join(files, "Alanin.xyz"))
    mol2 = Structure(os.path.join(files, "Alanin_different_atom_order_methyl_rotated_60_deg.xyz"))
    analyze = Analyzer()
    assert analyze.doubles(mol1, mol2, loose=True)


def test_strict_doubles_check():
    mol1 = Structure(os.path.join(files, "Tyrosin.xyz"))
    mol2 = Structure(os.path.join(files, "Tyrosin_double.xyz"))
    mol3 = Structure(os.path.join(files, "Tyrosin_phenyl_rotated_180_deg.xyz"))
    analyze = Analyzer()
    mol1_mol2 = analyze.doubles(mol1, mol2)
    mol1_mol3 = analyze.doubles(mol1, mol3)
    mol2_mol3 = analyze.doubles(mol2, mol3)
    assert (mol1_mol2 and mol1_mol3 and mol2_mol3)
