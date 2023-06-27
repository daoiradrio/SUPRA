#!/usr/bin/supra-python3



import os
import argparse
from SUPRAConformer.structure import Structure
from SUPRAConformer.conformergenerator import ConformerGenerator
from utils.analyzer import Analyzer



def main():
    mol = Structure()
    generator = ConformerGenerator()
    analyzer = Analyzer()
    parser = argparse.ArgumentParser()

    parser.add_argument("-path", type=str, required=True)
    args = parser.parse_args()

    mol.get_structure(os.path.abspath(args.path))

    print()

    generator.generate_conformers(mol)

    print("Performing removal of duplicate structures...")
    conformers = analyzer.remove_doubles(path=os.path.join(os.getcwd(), "SUPRA_Output/"), use_energy=True, mode="normal")
    conformers = analyzer.remove_doubles(path=os.path.join(os.getcwd(), "SUPRA_Output/"), use_energy=True, mode="tight")
    print("Removal of double structures done.")
    print(f"Individual generated conformers: {conformers}")

    print()


if __name__ == "__main__":
    main()
