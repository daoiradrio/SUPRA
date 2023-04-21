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
    A1 = Structure(os.path.join(files, "Alanin.xyz"))
    A2 = Structure(os.path.join(files, "Alanin_different_atom_order.xyz"))
    A3 = Structure(os.path.join(files, "Alanin_different_atom_order_methyl_rotated_60_deg.xyz"))
    T1 = Structure(os.path.join(files, "Tyrosin.xyz"))
    T2 = Structure(os.path.join(files, "Tyrosin_double.xyz"))
    T3 = Structure(os.path.join(files, "Tyrosin_phenyl_rotated_180_deg.xyz"))
    analyze = Analyzer()
    A1_A2 = analyze.doubles(A1, A2)
    A1_A3 = analyze.doubles(A1, A3)
    A2_A3 = analyze.doubles(A2, A3)
    T1_T2 = analyze.doubles(T1, T2)
    T1_T3 = analyze.doubles(T1, T3)
    T2_T3 = analyze.doubles(T2, T3)
    assert (A1_A2 and not A1_A3 and not A2_A3 and T1_T2 and T1_T3 and T2_T3)
