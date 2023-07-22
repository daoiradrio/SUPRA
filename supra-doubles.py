#!/usr/bin/supra-python3



import os
import argparse
from utils.analyzer import Analyzer
from SUPRAConformer.structure import Structure



def main():
    parser = argparse.ArgumentParser()
    analyzer = Analyzer()

    parser.add_argument("-path", type=str, required=False)
    parser.add_argument("-rmsd", type=float, required=False, default=0.1)
    parser.add_argument("-ignore", type=str, required=False, choices=["terminal", "methyl"], default=None)
    parser.add_argument("-path1", type=str, required=False)
    parser.add_argument("-path2", type=str, required=False)
    parser.add_argument("-matching", type=str, required=False, choices=["loose", "normal", "tight"], default="normal")
    args = parser.parse_args()
    
    if args.path:
        print()
        conformers_before = analyzer._count_conformers(path=args.path)
        print(f"Number of conformers in {args.path}: {conformers_before}")
        print("Performing removal of duplicate structures...")
        conformers_after = analyzer.remove_doubles(path=args.path, rmsd_threshold=args.rmsd, ignore=args.ignore, matching=args.matching)
        print("Removal of double structures done.")
        print(f"Individual conformers in {args.path}: {conformers_after}")
        print()
    else:
        if os.path.isfile(args.path1) and os.path.isfile(args.path2):
            structure1 = Structure(args.path1)
            structure2 = Structure(args.path2)
            if args.matching == "loose":
                kabsch_coords1, kabsch_coords2 = analyzer.kabsch(structure1.coords, structure2.coords)
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
            print("Comparing structures...")
            path1, n_conformers1, path2, n_conformers2, overlap = analyzer.compare_structure_sets(
                path1=args.path1,
                path2=args.path2,
                rmsd_threshold=args.rmsd,
                ignore=args.ignore,
                matching=args.matching
            )
            print("Comparing structures done.")
            print()
            print(f"Path 1: {path1}")
            print(f"Path 2: {path2}")
            print()
            print(f"Number of structures in Path 1: {n_conformers1}")
            print(f"Number of structures in Path 2: {n_conformers2}")
            print(f"Number of structures of Path 1 in Path 2: {overlap}")
            print()



if __name__ == "__main__":
    main()
