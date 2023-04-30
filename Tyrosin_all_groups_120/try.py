from SUPRAConformer.structure import Structure
from SUPRAConformer.conformergenerator import ConformerGenerator
import time

filename = "/home/baum/SUPRA/inputfiles/Tyrosin.xyz"
mol = Structure(filename)
gen = ConformerGenerator()
start = time.time()
gen.generate_conformers(mol)
print(f"Total time: {time.time() - start}")
