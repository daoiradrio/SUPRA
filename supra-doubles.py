#!/usr/bin/supra-python3



import os
import argparse
from utils.analyzer import Analyzer



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
        if os.path.isfile(args.path):
            conformers_before = analyzer._count_conformers_file(path=args.path)
            print()
            print(f"Number of structures in {args.path}: {conformers_before}")
            print("Performing removal of duplicate structures...")
            conformers_after = analyzer.remove_doubles_file(
                ensemble_file=args.path,
                rmsd_threshold=args.rmsd,
                ignore=args.ignore,
                matching=args.matching
            )
        else:
            conformers_before = analyzer._count_conformers_dir(path=args.path)
            print()
            print(f"Number of structures in {args.path}: {conformers_before}")
            print("Performing removal of duplicate structures...")
            conformers_after = analyzer.remove_doubles_dir(path=args.path, rmsd_threshold=args.rmsd, ignore=args.ignore, matching=args.matching)
        print("Removal of double structures done.")
        print(f"Individual conformers in {args.path}: {conformers_after}")
        print()
    elif args.path1 and args.path2:
        print()
        print("Comparing structures...")
        if os.path.isfile(args.path1) and os.path.isfile(args.path2):
            path1, n_conformers1, path2, n_conformers2, overlap = analyzer.compare_ensembles_files(
                path1=args.path1,
                path2=args.path2,
                rmsd_threshold=args.rmsd,
                ignore=args.ignore,
                matching=args.matching
            )
        else:
            path1, n_conformers1, path2, n_conformers2, overlap = analyzer.compare_ensembles_dirs(
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
    else:
        print()
        print("INVALID INPUT PATH(S)")
        print("SUPRA ABORTING.")
        print()



if __name__ == "__main__":
    main()
