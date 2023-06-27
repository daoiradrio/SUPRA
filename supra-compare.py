#!/usr/bin/supra-python3



import os
import argparse
from utils.analyzer import Analyzer



def main():
    parser = argparse.ArgumentParser()
    analyzer = Analyzer()

    parser.add_argument("-path1", type=str, required=True)
    parser.add_argument("-path2", type=str, required=True)
    parser.add_argument("-rmsd", type=float, required=False, const=0.1)
    parser.add_argument("-ignore", type=str, required=False, const=None)
    parser.add_argument("-mode", type=str, required=False, const="normal")
    args = parser.parse_args()

    analyzer.compare_structure_sets(path1=args.path1, path2=args.path2, rmsd_threshold=args.rmsd, ignore=args.ignore, mode=args.mode)



if __name__ == "__main__":
    main()
