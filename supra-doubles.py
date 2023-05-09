#!/usr/bin/supra-python3



import os
from utils.analyzer import Analyzer



def main():
    an = Analyzer()

    if len(os.sys.argv) == 3:
        if os.sys.argv[2].isnumeric():
            rmsd_threshold = float(os.sys.argv[2])
            strucs_path = os.sys.argv[1]
            os.path.abspath(strucs_path)
            print()
            an.remove_doubles(strucs_path, rmsd_threshold)
            print()
        else:
            print()
            print("USAGE: supra-doubles /path/to/structures/directory")
            print("OPTIONAL: supra-doubles /path/to/structures/directory [RMSD-Threshold]")
            print()
            print("[RMSD-Threshold] must be a number!")
            print()
            return
    elif len(os.sys.argv) == 2:
        strucs_path = os.sys.argv[1]
        os.path.abspath(strucs_path)
        print()
        an.remove_doubles(strucs_path)
        print()
    else:
        print()
        print("USAGE: supra-doubles /path/to/structures/directory")
        print("OPTIONAL: supra-doubles /path/to/structures/directory [RMSD-Threshold]")
        print()
        return



if __name__ == "__main__":
    main()
