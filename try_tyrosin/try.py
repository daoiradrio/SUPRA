from SUPRAConformer.structure import Structure
from SUPRAConformer.conformergenerator import ConformerGenerator

filename = "/home/baum/SUPRA/inputfiles/Tyrosin.xyz"
mol = Structure(filename)
gen = ConformerGenerator()
gen.generate_conformers(mol)
