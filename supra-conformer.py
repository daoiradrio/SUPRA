#!/usr/bin/supra-python3



import os
from SUPRAConformer.structure import Structure
from SUPRAConformer.conformergenerator import ConformerGenerator
from utils.analyzer import Analyzer



def main():
    if len(os.sys.argv) != 2:
        print("\nUSAGE: supra-conformer /path/to/file\n")
        return
    else:
        filepath = os.sys.argv[1]

    mol = Structure()
    gen = ConformerGenerator()
    an = Analyzer()

    print()
    mol.get_structure(filepath)
    gen.generate_conformers(mol)
    an.remove_doubles("SUPRA_Output/")
    print()


if __name__ == "__main__":
    main()
