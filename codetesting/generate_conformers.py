from SUPRAConformer.structure import Structure
from SUPRAConformer.conformergenerator import ConformerGenerator



filepath = "inputfiles/Tyrosin.xyz"
mol = Structure(filepath)
generator = ConformerGenerator()
generator.generate_conformers(mol)
