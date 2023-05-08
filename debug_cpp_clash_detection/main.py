from SUPRAConformer.structure import Structure
from SUPRAConformer.conformergenerator import ConformerGenerator
import os

filepath = "/home/baum/SUPRA/inputfiles/Tyrosin.xyz"
#clash_filepath = os.path.join("debug_files", f"clash_structure{os.sys.argv[1]}.xyz")
clash_filepath = os.path.join("debug_files", "clash_structure1.xyz")

mol = Structure()
clash_mol = Structure()
gen = ConformerGenerator()

mol.get_structure(filepath)
clash_mol.get_structure(clash_filepath)
if gen._clashes(mol.bond_partners, clash_mol.coords):
    print("CLASHES!")
else:
    print("NO CLASHES.")

#print("___")
#n = len(os.listdir("debug_files"))
#for i in range(n):
#    clash_filepath = os.path.join("debug_files", f"clash_structure{i}.xyz")
#    clash_mol.get_structure(clash_filepath)
#    if gen._clashes(mol.bond_partners, clash_mol.coords):
#        print("CLASHES!")
#    else:
#        print("NO CLASHES.")
#print("___")
