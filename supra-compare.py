#!/usr/bin/supra-python3



import os
import argparse
from utils.analyzer import Analyzer



def main():
    parser = argparse.ArgumentParser()
    analyzer = Analyzer()

    parser.add_argument("-path1", type=str, required=True)
    parser.add_argument("-path2", type=str, required=True)
    parser.add_argument("-rmsd", type=float, required=False)
    parser.add_argument("-ignore", type=str, required=False)
    args = parser.parse_args()

    if args.rmsd:
        rmsd = args.rmsd
    else:
        rmsd = 0.1

    analyzer.compare_structure_sets(path1=args.path1, path2=args.path2, rmsd_threshold=rmsd, ignore=args.ignore)



if __name__ == "__main__":
    main()
