#!/usr/bin/supra-python3



import os
import argparse
from SUPRAConformer.structure import Structure
from SUPRAConformer.conformergenerator import ConformerGenerator
from utils.analyzer import Analyzer



def main():
    mol = Structure()
    gen = ConformerGenerator()
    an = Analyzer()
    parser = argparse.ArgumentParser()

    parser.add_argument("-path", type=str, required=True)
    args = parser.parse_args()

    print()
    mol.get_structure(os.path.abspath(args.path))
    gen.generate_conformers(mol)
    an.remove_doubles(path=os.path.join(os.getcwd(), "SUPRA_Output/"), use_energy=True)
    print()


if __name__ == "__main__":
    main()
