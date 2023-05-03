from utils.analyzer import Analyzer
from SUPRAConformer.structure import Structure
from SUPRAConformer.conformergenerator import ConformerGenerator
import time
import os

supra_folder = "/home/baum/SUPRA"
file = "inputfiles/Tryptophan.xyz"
output = "Tryptophan_all_groups_90/SUPRA_Output"
filename = os.path.join(supra_folder, file)
output_folder = os.path.join(supra_folder, output)

mol = Structure()
gen = ConformerGenerator()

mol.get_structure(filename)
start = time.time()
gen.generate_conformers(mol)
time_generation = time.time() - start
start = time.time()
Analyzer.remove_doubles(output_folder)
time_analyzation = time.time() - start
print(f"Time generation: {time_generation}")
print(f"Time analyzation: {time_analyzation}")
