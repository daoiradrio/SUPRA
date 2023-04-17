import os
import shutil
import subprocess
import importlib.util
#from importlib.util import find_spec

import numpy as np

from Helper import generate_folder, guess_molecule_name
from Converter import convert_xyz_to_smiles
from rdkit import Chem
from rdkit.Chem import rdMolTransforms
from scipy import ndimage
from PeriodicTable import atomic_masses_symbols



class Geometry:
        
    def __init__(self):
        pass
    

    @staticmethod
    def distance(pos1, pos2):
        vec1 = np.array(pos1)
        vec2 = np.array(pos2)
        return np.linalg.norm(vec1 - vec2)



class ParameterGenerator:
    
    def __init__(self, infile, folder):
        self.infile = os.path.abspath(infile)
        self.folder = folder
        self.molecule_name = None
        self.smiles = None
        self.xyz_coordinates = None

        abs_path_new_folder = os.path.abspath(folder)
        if os.getcwd != abs_path_new_folder:
            if os.path.exists(abs_path_new_folder):
                shutil.rmtree(abs_path_new_folder, ignore_errors=True)
            os.makedirs(abs_path_new_folder)

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
            print("\nERROR! BOSS is not found on your work station. $BOSSdir path is not set. DO NOT forget to export the BOSSdir path (export BOSSdir=/path/boss/intalled)\n")
            exit()
        if importlib.util.find_spec('rdkit') == None:
            print("\nERROR! RDKit is not installed on your work station.\n")
            exit()
        
        self.smiles, self.xyz_coordinates = convert_xyz_to_smiles(self.infile)
        mol = Chem.MolFromSmiles(self.smiles)
        properly_ordered_atoms_indices_list = self.properly_ordered_atoms_list(mol, abs_path_new_folder)
        mol, newIndexToOriginalIndex = self.mol_with_proper_atom_order(mol, properly_ordered_atoms_indices_list)

    
    def properly_ordered_atoms_list(self, mol, workdir):
        heavy_atoms_indices = [atom.GetIdx() for atom in mol.GetAtoms() if atom.GetSymbol() != 'H']
        hydrogen_atoms_indices = [atom.GetIdx() for atom in mol.GetAtoms() if atom.GetSymbol() == "H"]
        index_atom_closest_to_com = self.atom_index_closest_to_com(mol)
        atoms_indices_with_proper_order = [index_atom_closest_to_com]
        self.add_atoms_to_index_list(mol, atoms_indices_with_proper_order, heavy_atoms_indices)
        self.add_atoms_to_index_list(mol, atoms_indices_with_proper_order, hydrogen_atoms_indices)
        return atoms_indices_with_proper_order


    def atom_index_closest_to_com(self, mol):
        #***
        masses = [atomic_masses_symbols[atom.GetSymbol()] for atom in mol.GetAtoms()]
        coords_and_masses = [np.array([coords[0], coords[1], coords[2], mass]) for coords, mass in zip(self.xyz_coordinates, masses)]
        coords_and_masses = np.array(coords_and_masses)
        com = np.average(coords_and_masses, axis=0, weights=coords_and_masses[:,3])
        com = com[:3]
        # ***
        distances = []
        for atom, atom_position in zip(mol.GetAtoms(), self.xyz_coordinates):
            if atom.GetSymbol() == 'H':
                continue
            distances.append([Geometry.distance(atom_position, com), atom])
        atom_closest_to_com = sorted(distances, key=lambda x: x[0])[0]
        return atom_closest_to_com[1].GetIdx()


    def add_atoms_to_index_list(self, mol, atoms_indices_with_proper_order, heavy_atoms_indices):
        missing_atoms = []
        for heavy_atom_index in heavy_atoms_indices:
            heavy_atom = mol.GetAtomWithIdx(heavy_atom_index)
            neighbors_indices = [x.GetIdx() for x in heavy_atom.GetNeighbors()]
            if heavy_atom_index not in atoms_indices_with_proper_order:
                if any([True for neighbor_index in neighbors_indices if neighbor_index in atoms_indices_with_proper_order]):
                    atoms_indices_with_proper_order.append(heavy_atom.GetIdx())
                else:
                    missing_atoms.append(heavy_atom.GetIdx())
        if len(missing_atoms) == 0:
            return
        else:
            self.add_atoms_to_index_list(mol, atoms_indices_with_proper_order, heavy_atoms_indices)
        

    def mol_with_proper_atom_order(self, mol, atoms_indices_list):
        new_atoms_labels = {}
        new_to_original_index = {}
        new_molecule_block = ""
        for i, index in enumerate(atoms_indices_list):
            new_to_original_index[i] = index
            atom = mol.GetAtomWithIdx(index)
            atom_position = self.xyz_coordinates[index]
            try:
                print("Hier")
                residue_info = atom.GetPDBResidueInfo()
                residue_name = residue_info.GetResidueName()
                atom_name = residue_info.GetName()
            except:
                print("Hmmm...")
                element = atom.GetSymbol()
                if not element in new_atoms_labels.keys():
                    new_atoms_labels[element] = 1
                else:
                    new_atoms_labels[element] += 1
                atom_name = f"{element}{new_atoms_labels[element]}"
                residue_name = "MOL"
            line = 'ATOM%7d%5s%4s%6d%12.3f%8.3f%8.3f  1.00  0.00%12s  \n' % (i, atom_name, residue_name, 1, atom_position[0], atom_position[1], atom_position[2], atom.GetSymbol())
            new_molecule_block += line
        print(new_molecule_block)
        new_mol = Chem.MolFromPDBBlock(new_molecule_block, removeHs=False)
        return new_mol, new_to_original_index



input_xyz = os.sys.argv[1]
workdir = "test"
ParameterGenerator(input_xyz, workdir)
