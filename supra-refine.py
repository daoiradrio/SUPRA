#!/usr/bin/supra-python3



import os
from SUPRAConformer.optimizer import Optimizer
from utils.analyzer import Analyzer



def main():
    if len(os.sys.argv) != 2:
        print("\nUSAGE: supra-refine /path/to/structures/directory\n")
        return
    else:
        strucs_path = os.sys.argv[1]
        os.path.abspath(strucs_path)

    opt = Optimizer()
    an = Analyzer()

    opt.refine_ff_opts(strucs_path)
    an.remove_doubles(strucs_path)



if __name__ == "__main__":
    main()
