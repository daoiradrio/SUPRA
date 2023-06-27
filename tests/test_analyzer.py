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
    rmsd_threshold = 0.1
    analyze = Analyzer()
    mol1 = Structure(os.path.join(files, "Alanin.xyz"))
    mol2 = Structure(os.path.join(files, "Alanin_methyl_rotated_60_deg.xyz"))

    for atom in analyze.get_methyl_group_atoms(mol1.bond_partners):
        del mol1.coords[atom]
    for atom in analyze.get_methyl_group_atoms(mol2.bond_partners):
        del mol2.coords[atom]

    assert analyze.rmsd(mol1.coords, mol2.coords) <= rmsd_threshold



def test_strict_doubles_check():
    rmsd_threshold = 0.1
    analyze = Analyzer()
    #A1 = Structure(os.path.join(files, "Alanin.xyz"))
    #A2 = Structure(os.path.join(files, "Alanin_different_atom_order.xyz"))
    #A3 = Structure(os.path.join(files, "Alanin_different_atom_order_methyl_rotated_60_deg.xyz"))
    T1 = Structure(os.path.join(files, "Tyrosin.xyz"))
    T2 = Structure(os.path.join(files, "Tyrosin_double.xyz"))
    T3 = Structure(os.path.join(files, "Tyrosin_phenyl_rotated_180_deg.xyz"))
    T4 = Structure(os.path.join(files, "Alanin.xyz"))
    T5 = Structure(os.path.join(files, "Alanin_methyl_rotated_60_deg.xyz"))

    #A1_A2 = analyze.doubles(A1.coords, A2.coords, rmsd_threshold)
    #A1_A3 = analyze.doubles(A1.coords, A3.coords, rmsd_threshold)
    #A2_A3 = analyze.doubles(A2.coords, A3.coords, rmsd_threshold)
    T1_T2 = analyze.rmsd(T1.coords, T2.coords) <= rmsd_threshold
    T1_T3 = analyze.rmsd(T1.coords, T3.coords) <= rmsd_threshold
    T2_T3 = analyze.rmsd(T2.coords, T3.coords) <= rmsd_threshold
    T4_T5 = analyze.rmsd(T4.coords, T5.coords) <= rmsd_threshold
    
    #assert (A1_A2 and not A1_A3 and not A2_A3 and T1_T2 and T1_T3 and T2_T3)
    assert (T1_T2 and T1_T3 and T2_T3 and not T4_T5)



def test_compare_structure_sets():
    rmsd_threshold = 0.1
    analyze = Analyzer()
    set1 = os.path.join(files, "Alanin_set1")
    set2 = os.path.join(files, "Alanin_set2")

    strict = analyze.compare_structure_sets(
                path1=set1,
                path2=set2,
                rmsd_threshold=rmsd_threshold
             )

    loose = analyze.compare_structure_sets(
                path1=set1,
                path2=set2,
                rmsd_threshold=rmsd_threshold,
                ignore="methyl"
            )

    assert (strict == 1 and loose == 2)
