from SUPRAConformer.structure import Structure
from utils.analyzer import Analyzer



struc1 = Structure("inputfiles/Alanin.xyz")
struc2 = Structure("inputfiles/Alanin_rotated_methyl.xyz")
an = Analyzer()

an.kabsch(struc1.coords, struc2.coords)
print(an.kabsch_and_rmsd(struc1.coords, struc2.coords))
