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
    parser.add_argument("-chrg", type=int, required=False)
    args = parser.parse_args()
    strucs_path = os.path.abspath(args.path)
    if args.chrg:
        chrg = args.chrg
    else:
        chrg = None

    #print()
    #opt.refine_structures_xtb(strucs_path, chrg)
    #analyzer.remove_doubles(strucs_path)
    #print()

    print("Performing removal of duplicate structures...")
    conformers = analyzer.remove_doubles(path=os.path.join(os.getcwd(), "SUPRA_Output/"), use_energy=True, mode="normal")
    conformers = analyzer.remove_doubles(path=os.path.join(os.getcwd(), "SUPRA_Output/"), use_energy=True, mode="tight")
    print("Removal of double structures done.")
    print(f"Individual generated conformers: {conformers}")



if __name__ == "__main__":
    main()
