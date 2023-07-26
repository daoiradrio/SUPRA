#!/usr/bin/supra-python3



import os
import argparse
from utils.analyzer import Analyzer
from SUPRAConformer.structure import Structure



def main():
    parser = argparse.ArgumentParser()
    analyzer = Analyzer()

    parser.add_argument("-path1", type=str, required=True)
    parser.add_argument("-path2", type=str, required=True)
    parser.add_argument("-ignore", type=str, required=False, choices=["terminal", "methyl"], default=None)
    parser.add_argument("-matching", type=str, required=False, choices=["loose", "normal", "tight"], default="normal")
    args = parser.parse_args()


    if os.path.isfile(args.path1) and os.path.isfile(args.path2):
        structure1 = Structure(args.path1)
        structure2 = Structure(args.path2)
        if args.matching == "loose":
            kabsch_coords1, kabsch_coords2 = analyzer._kabsch(structure1.coords, structure2.coords)
            rmsd = analyzer._calc_rmsd(kabsch_coords1, kabsch_coords2)
        elif args.matching == "normal":
            rmsd = analyzer._rmsd(structure1.coords, structure2.coords)
        elif args.matching == "tight":
            rmsd = analyzer._rmsd_tight(structure1.coords, structure1.bond_partners, structure2.coords, structure2.bond_partners)
        print()
        print(f"Path of molecule 1: {args.path1}")
        print(f"Path of molecule 2: {args.path2}")
        print()
        print(f"RMSD: {rmsd:.4f}")
        print()
    else:
        print()
        print("INVALID INPUT FILE(S):")
        print(f"Input path 1: {args.path1}")
        print(f"Input path 2: {args.path2}")
        print("SUPRA ABORTING.")
        print()



if __name__ == "__main__":
    main()
