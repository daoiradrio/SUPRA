#!/usr/bin/supra-python3



import os
import argparse
from utils.optimizer import Optimizer
from utils.analyzer import Analyzer



def main():
    opt = Optimizer()
    analyzer = Analyzer()
    parser = argparse.ArgumentParser()

    parser.add_argument("-path", type=str, required=True)
    parser.add_argument("-code", type=str, required=False, choices=["pyscf", "xtb"], default="xtb")
    parser.add_argument("-chrg", type=int, required=False, default=None)
    args = parser.parse_args()

    print()

    print("Performining refining structure optimization...")
    if args.code == "pyscf":
        opt.qc_refine_structures(os.path.abspath(args.path))
    elif args.code == "xtb":
        opt.xtb_refine_structures(os.path.abspath(args.path), args.chrg)
    print("Refining structure optimization done.")

    print("Performing removal of duplicate structures...")
    conformers = analyzer.remove_doubles(path=args.path, use_energy=False, matching="normal")
    conformers = analyzer.remove_doubles(path=args.path, use_energy=False, matching="tight")
    print("Removal of double structures done.")
    print(f"Individual generated conformers: {conformers}")
    
    print()



if __name__ == "__main__":
    main()
