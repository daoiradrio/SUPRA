import os
from utils.analyzer import Analyzer

curr_dir = os.getcwd()
workdir = os.path.join(curr_dir, "conformers")
struc_filename = "conformer"

os.makedirs(workdir)

with open("Alanin_ensemble.xyz", "r") as infile:
    n_atoms = int(infile.readline().split()[0])
    file_counter = 0
    line_iter = 1
    new_struc = os.path.join(workdir, f"{struc_filename}{file_counter}.xyz")
    new_struc_file = open(new_struc, "w")
    print(n_atoms, file=new_struc_file)
    for line in infile:
        if (line_iter < n_atoms+2):
            print(line, file=new_struc_file, end="")
            line_iter += 1
        else:
            new_struc_file.close()
            file_counter += 1
            line_iter = 1
            new_struc = os.path.join(workdir, f"{struc_filename}{file_counter}.xyz")
            new_struc_file = open(new_struc, "w")
            print(n_atoms, file=new_struc_file)

analyzer = Analyzer()
analyzer.remove_doubles(path=workdir, rmsd_threshold=0.1, ignore="methyl")