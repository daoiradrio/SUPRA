#!/usr/bin/supra-python3



import os
import argparse
from utils.analyzer import Analyzer



def main():
    parser = argparse.ArgumentParser()
    an = Analyzer()

    parser.add_argument("-path", type=str, required=True)
    parser.add_argument("-rmsd", type=float, required=False)
    parser.add_argument("-ignore", type=str, required=False)
    args = parser.parse_args()

    print()
    an.remove_doubles(path=os.path.abspath(args.path), rmsd_threshold=args.rmsd, ignore=args.ignore)
    print()



if __name__ == "__main__":
    main()
