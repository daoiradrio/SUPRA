#!/usr/bin/supra-python3



import os
from utils.analyzer import Analyzer



# TODO: UPDATE WITH ARGPARSE PACKAGE!
def main():
    an = Analyzer()
    if len(os.sys.argv) == 4:
        if os.sys.argv[3] == "loose":
            try:
                rmsd_threshold = float(os.sys.argv[2])
                strucs_path = os.sys.argv[1]
                os.path.abspath(strucs_path)
                print()
                an.remove_doubles(path=strucs_path, rmsd_threshold=rmsd_threshold, loose=True)
                print()
            except:
                print()
                print("USAGE: supra-doubles /path/to/structures/directory")
                print("OPTIONAL: supra-doubles /path/to/structures/directory [RMSD-Threshold]")
                print("OPTIONAL: supra-doubles /path/to/structures/directory [RMSD-Threshold] loose")
                print()
                return
        else:
            print()
            print("USAGE: supra-doubles /path/to/structures/directory")
            print("OPTIONAL: supra-doubles /path/to/structures/directory [RMSD-Threshold]")
            print("OPTIONAL: supra-doubles /path/to/structures/directory [RMSD-Threshold] loose")
            print()
            return
    elif len(os.sys.argv) == 3:
        if os.sys.argv[2] == "loose":
            strucs_path = os.sys.argv[1]
            os.path.abspath(strucs_path)
            print()
            an.remove_doubles(path=strucs_path, loose=True)
            print()
        else:
            try:
                rmsd_threshold = float(os.sys.argv[2])
                strucs_path = os.sys.argv[1]
                os.path.abspath(strucs_path)
                print()
                an.remove_doubles(path=strucs_path, rmsd_threshold=rmsd_threshold)
                print()
            except:
                print()
                print("USAGE: supra-doubles /path/to/structures/directory")
                print("OPTIONAL: supra-doubles /path/to/structures/directory [RMSD-Threshold]")
                print("OPTIONAL: supra-doubles /path/to/structures/directory [RMSD-Threshold] loose")
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
