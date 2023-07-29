import os

from SUPRAConformer.structure import Structure
from utils.analyzer import Analyzer



files = "testcases/"
try:
    f = open(os.path.join(files, "Alanin.xyz"))
    f.close()
except:
    files = "tests/testcases/"



def test_doubles_check_ignore_methyl():
    rmsd_threshold = 0.1
    analyze = Analyzer()
    mol1 = Structure(os.path.join(files, "Alanin.xyz"))
    mol2 = Structure(os.path.join(files, "Alanin_methyl_rotated_60_deg.xyz"))

    for atom in analyze._find_methyl_group_atoms(mol1.bond_partners):
        del mol1.coords[atom]
    for atom in analyze._find_methyl_group_atoms(mol2.bond_partners):
        del mol2.coords[atom]

    assert analyze._rmsd(mol1.coords, mol2.coords) <= rmsd_threshold



def test_doubles_check_hungarian_atom_matching():
    rmsd_threshold = 0.1
    analyze = Analyzer()

    T1 = Structure(os.path.join(files, "Tyrosin.xyz"))
    T2 = Structure(os.path.join(files, "Tyrosin_double.xyz"))
    T3 = Structure(os.path.join(files, "Tyrosin_phenyl_rotated_180_deg.xyz"))
    T4 = Structure(os.path.join(files, "Alanin.xyz"))
    T5 = Structure(os.path.join(files, "Alanin_methyl_rotated_60_deg.xyz"))

    T1_T2 = analyze._rmsd(T1.coords, T2.coords) <= rmsd_threshold
    T1_T3 = analyze._rmsd(T1.coords, T3.coords) <= rmsd_threshold
    T2_T3 = analyze._rmsd(T2.coords, T3.coords) <= rmsd_threshold
    T4_T5 = analyze._rmsd(T4.coords, T5.coords) <= rmsd_threshold
    
    assert (T1_T2 and T1_T3 and T2_T3 and not T4_T5)



def test_doubles_check_tight_matching():
    rmsd_threshold = 0.1
    analyzer = Analyzer()

    A1 = Structure(os.path.join(files, "Alanin.xyz"))
    A2 = Structure(os.path.join(files, "Alanin_different_atom_order.xyz"))
    P1 = Structure(os.path.join(files, "Tyr-Ala-Trp.xyz"))
    P2 = Structure(os.path.join(files, "Tyr-Ala-Trp_different_order.xyz"))

    A1_A2 = analyzer._rmsd_tight(A1.coords, A1.bond_partners, A2.coords, A2.bond_partners) <= rmsd_threshold
    P1_P2 = analyzer._rmsd_tight(P1.coords, P1.bond_partners, P2.coords, P2.bond_partners) <= rmsd_threshold

    assert A1_A2



def test_compare_structure_ensemble_dirs():
    rmsd_threshold = 0.1
    analyze = Analyzer()
    set1 = os.path.join(files, "Alanin_set1")
    set2 = os.path.join(files, "Alanin_set2")

    _, _, _, _, strict = analyze.compare_ensembles_dirs(
                path1=set1,
                path2=set2,
                rmsd_threshold=rmsd_threshold
             )

    _, _, _, _, loose = analyze.compare_ensembles_dirs(
                path1=set1,
                path2=set2,
                rmsd_threshold=rmsd_threshold,
                ignore="methyl"
            )

    assert (strict == 1 and loose == 2)



def test_check_for_duplicates():
    rmsd_threshold = 0.1
    analyzer = Analyzer()

    xyz_file = os.path.join(files, "Alanin.xyz")
    set1 = os.path.join(files, "Alanin_set1")
    set2 = os.path.join(files, "Alanin_set2")

    res1 = analyzer.check_for_duplicates(xyz_file=xyz_file, path=set1, matching="tight")
    res2 = analyzer.check_for_duplicates(xyz_file=xyz_file, path=set2, matching="tight")

    assert (res1 == xyz_file and res2 == None)