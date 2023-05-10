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
    args = parser.parse_args()
    strucs_path = os.path.abspath(args.path)

    print()
    opt.refine_ff_opts(strucs_path)
    an.remove_doubles(strucs_path)
    print()



if __name__ == "__main__":
    main()
