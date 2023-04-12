import os
import shutil
import subprocess
import importlib.util

import numpy as np

from Helper import generate_folder, guess_molecule_name
from rdkit import Chem
from rdkit.Chem import rdMolTransforms



class Geometry:
        
    def __init__(self):
        pass
    

    @staticmethod
    def distance(pos1: list, pos2: list):
        vec1 = np.array(pos1)
        vec2 = np.array(pos2)
        return np.linalg.norm(vec1 - vec2)



class ParameterGenerator:
    
    def __init__(self, infile, folder, molecule_name = None):
        self.infile = os.path.abspath(infile)
        self.folder = folder
        self.molecule_name = molecule_name

        abs_path_new_folder = os.path.abspath(folder)
        if os.getcwd != abs_path_new_folder:
            if os.path.exists(abs_path_new_folder):
                shutil.rmtree(abs_path_new_folder, ignore_errors=True)
            os.makedirs(abs_path_new_folder)

        if molecule_name == None:
            molecule_name = os.path.splitext(os.path.basename(infile))[0]

        babel = shutil.which("obabel")
        if babel == None:
            print("\nERROR! Open Babel is not installed on your work station.\n")
            exit()
        csh = shutil.which("csh")
        if csh == None:
            print("\nERROR! csh is not installed on your work station.\n")
            exit()
        try:
            os.environ["BOSSdir"]
        except:
            print("\nERROR! BOSS is not found on your work station. $BOSSdir path is not set. DO NOT forget to export the BOSSdir path (export BOSSdir=/path/boss/intalled)")
            exit()
        if importlib.util.find_spec('rdkit') == None:
            print("\nERROR! RDKit is not installed on your work station.\n")
            exit()

        mol = Chem.MolFromXYZFile(infile)
        properly_ordered_atoms_indices_list = self.properly_ordered_atoms_list(mol, molecule_name, abs_path_new_folder)
        #mol, newIndexToOriginalIndex = self.mol_with_proper_atom_order(mol, properly_ordered_atoms_indices_list)

    
    def properly_ordered_atoms_list(self, mol, molecule_name, workdir):
         heavy_atoms_indices = [atom.GetIdx() for atom in mol.GetAtoms() if atom.GetSymbol() != 'H'] 
         index_atom_closest_to_com = self.atom_index_closest_to_com(mol)


    def atom_index_closest_to_com(self, mol):
        conf = mol.GetConformer(0)
        com = rdMolTransforms.ComputeCentroid(conf)
        distances = []
        for atom in mol.GetAtoms():
            if atom.GetSymbol() == 'H':
                continue
            atom_position = conf.GetAtomPosition(atom.GetIdx())
            distances.append([Geometry.distance([atom_position.x, atom_position.y, atom_position.z], com), atom])
        atom_closest_to_com = sorted(distances, key=lambda x: x[0])[0]
        return atom_closest_to_com[1].GetIdx()


    def mol_with_proper_atom_order(self, mol, atoms_indices_list):
        pass



ParameterGenerator("Alanin.xyz", "test")
