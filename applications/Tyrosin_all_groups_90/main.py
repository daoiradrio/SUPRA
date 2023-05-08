from utils.analyzer import Analyzer
from SUPRAConformer.structure import Structure
from SUPRAConformer.conformergenerator import ConformerGenerator
import os

supra_folder = "/home/baum/SUPRA"
file = "inputfiles/Tyrosin.xyz"
output = "applications/Tyrosin_all_groups_90/SUPRA_Output"
filename = os.path.join(supra_folder, file)
output_folder = os.path.join(supra_folder, output)

mol = Structure()
gen = ConformerGenerator()
an = Analyzer()

mol.get_structure(filename)
gen.generate_conformers(mol)
an.remove_doubles(output_folder)
