#!/usr/bin/supra-python3



import os
import argparse
from utils.analyzer import Analyzer



def main():
    parser = argparse.ArgumentParser()
    an = Analyzer()

    parser.add_argument("-path1", type=str, required=True)
    parser.add_argument("-path2", type=str, required=True)
    parser.add_argument("-rmsd", type=float, required=False)
    parser.add_argument("-ignore", type=str, required=False)
    args = parser.parse_args()

    print()
    an.compare_structure_sets(args.path1, args.path2)
    print()



if __name__ == "__main__":
    main()
