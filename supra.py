#!/home/dario/mambaforge/envs/supra/bin/python3



import os
from SUPRAConformer.structure import Structure
from SUPRAConformer.conformergenerator import ConformerGenerator
from utils.analyzer import Analyzer



def main():
    if len(os.sys.argv) != 2:
        print("USAGE: supra.py /path/to/file")
        return
    else:
        filepath = os.sys.argv[1]

    mol = Structure()
    gen = ConformerGenerator()
    an = Analyzer()

    mol.get_structure(filepath)
    gen.generate_conformers(mol)


if __name__ == "__main__":
    main()
