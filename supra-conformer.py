#!/usr/bin/supra-python3



import os
import argparse
from SUPRAConformer.structure import Structure
from SUPRAConformer.conformergenerator import ConformerGenerator

from time import time



def main():
    mol = Structure()
    generator = ConformerGenerator()
    parser = argparse.ArgumentParser()

    parser.add_argument("-path", type=str, required=True)
    parser.add_argument(
        "-ignore",
        type=str,
        required=False,
        choices=["methyl", "terminal", "peptide"],
        default=[],
        nargs="*"
    )
    parser.add_argument("-increment", type=int, required=False, choices=[60, 90, 120, 180], default=120)

    args = parser.parse_args()

    ignore_methyl = False
    ignore_terminal = False
    ignore_peptide = False
    if "methyl" in args.ignore:
        ignore_methyl = True
    if "terminal" in args.ignore:
        ignore_terminal = True
    if "peptide" in args.ignore:
        ignore_peptide = True

    mol.get_structure(os.path.abspath(args.path))

    start = time()
    n_conformers = generator.new_generate_conformers(
        structure=mol,
        increment=args.increment,
        ignore_methyl=ignore_methyl,
        ignore_terminal=ignore_terminal,
        ignore_peptide=ignore_peptide
    )
    stop = time()
    print(f"Benötigte Zeit: {stop-start}")



if __name__ == "__main__":
    main()
