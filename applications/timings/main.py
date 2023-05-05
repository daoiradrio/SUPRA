from utils.analyzer import Analyzer
from SUPRAConformer.structure import Structure
from SUPRAConformer.conformergenerator import ConformerGenerator
import time
import os

supra_folder = "/home/baum/SUPRA"
file = "inputfiles/Tryptophan.xyz"
output = "applications/timings/SUPRA_Output"
filename = os.path.join(supra_folder, file)
output_folder = os.path.join(supra_folder, output)

mol = Structure()
gen = ConformerGenerator()
an = Analyzer()

mol.get_structure(filename)
start = time.time()
gen.generate_conformers(mol)
time_generation = time.time() - start
start = time.time()
an.remove_doubles(output_folder)
time_analyzation = time.time() - start
print(f"Time generation: {time_generation}")
print(f"Time analyzation: {time_analyzation}")
print(f"Time for for: {an.time_for_for_for_for}")
print(f"Time read xyz: {an.time_read_xyz}")
print(f"Time Kabsch: {an.time_kabsch}")
print(f"Time build cost matrix: {an.time_build_cost_matrix}")
print(f"Time get_element: {an.time_get_element}")
print(f"Time if else elements match: {an.time_if_else_elements_match}")
print(f"Time call cost function: {an.time_call_cost_function}")
print(f"Time assign cost values: {an.time_assign_cost_values}")
print(f"Time Hungarian: {an.time_hungarian}")
print(f"Time RMSD: {an.time_rmsd}")
print(f"Time remove: {an.time_remove}")
