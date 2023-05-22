#!/usr/bin/supra-python3



import os
import argparse
from SUPRAConformer.optimizer import Optimizer
from utils.analyzer import Analyzer



def main():
    opt = Optimizer()
    an = Analyzer()
    parser = argparse.ArgumentParser()

    parser.add_argument("-path", type=str, required=True)
    parser.add_argument("-chrg", type=int, required=False)
    args = parser.parse_args()
    strucs_path = os.path.abspath(args.path)
    if args.chrg:
        chrg = args.chrg
    else:
        chrg = None

    print()
    opt.refine_ff_opts(strucs_path, chrg)
    an.remove_doubles(strucs_path)
    print()



if __name__ == "__main__":
    main()
