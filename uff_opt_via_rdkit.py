from utils.helper import get_element
from SUPRAConformer.structure import Structure

from rdkit import Chem
from rdkit.Chem import rdDetermineBonds
from rdkit.Chem.rdForceFieldHelpers import UFFOptimizeMolecule, UFFOptimizeMoleculeConfs



structure = Structure("/home/dario/SUPRA/inputfiles/Alanin.xyz")
atoms = [get_element(atom) for atom in structure.coords.keys()]
coords = [coord for coord in structure.coords.values()]
xyz_string = ""
xyz_string = xyz_string + f"{len(atoms)}\n\n"
for atom, (x,y,z) in zip(atoms, coords):
    xyz_string = xyz_string + f"{atom}\t{x}\t{y}\t{z}\n"

mol = Chem.MolFromXYZBlock(xyz_string)
mol = Chem.Mol(mol)
rdDetermineBonds.DetermineBonds(mol)

#"""
print()
for i, atom in enumerate(mol.GetAtoms()):
    pos = mol.GetConformer().GetAtomPosition(i)
    print(f"{atom.GetSymbol()}\t{pos.x}\t{pos.y}\t{pos.z}")
print()
#"""

res = UFFOptimizeMoleculeConfs(mol, maxIters=1000)
energy = res[0][1]
print(energy)

#"""
print()
for i, atom in enumerate(mol.GetAtoms()):
    pos = mol.GetConformer().GetAtomPosition(i)
    print(f"{atom.GetSymbol()}\t{pos.x}\t{pos.y}\t{pos.z}")
print()
#"""