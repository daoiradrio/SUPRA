#!/usr/bin/supra-python3



import os
import argparse
from utils.analyzer import Analyzer
from SUPRAConformer.structure import Structure



def main():
    parser = argparse.ArgumentParser()
    analyzer = Analyzer()

    parser.add_argument("-path", type=str, required=True)
    parser.add_argument("-rmsd", type=float, required=False)
    parser.add_argument("-ignore", type=str, required=False)
    parser.add_argument("-path1", type=str, required=False)
    parser.add_argument("-path2", type=str, required=False)
    parser.add_argument("-mode", type=str, required=False)
    args = parser.parse_args()

    if args.rmsd:
        rmsd = args.rmsd
    else:
        rmsd = 0.1
    
    if args.mode:
        mode = args.mode
    else:
        mode = "normal"
    
    if args.path:
        print()
        print("Performing removal of duplicate structures...")
        conformers = analyzer.remove_doubles(path=os.path.abspath(args.path), rmsd_threshold=rmsd, ignore=args.ignore)
        print("Removal of double structures done.")
        print(f"Individual conformers in {args.path}: {conformers}")
        print()
    else:
        structure1 = Structure(args.path1)
        structure2 = Structure(args.path2)
        print()
        if mode == "loose":
            rmsd = analyzer.calc_rmsd(structure1.coords, structure2.coords)
        elif mode == "normal":
            rmsd = analyzer.rmsd(structure1.coords, structure2.coords)
        elif mode == "tight":
            rmsd = analyzer.rmsd_tight(structure1.coords, structure1.bond_partners, structure2.coords, structure2.bond_partners)
        print(f"Path of molecule 1: {args.path1}")
        print(f"Path of molecule 2: {args.path2}")
        print()
        print(f"RMSD: {rmsd}")
        print()



if __name__ == "__main__":
    main()
