import os
import subprocess

from utils.helper import get_element
from rdkit import Chem
from rdkit.Chem import rdDetermineBonds
from rdkit.Chem.rdForceFieldHelpers import UFFOptimizeMoleculeConfs, MMFFOptimizeMoleculeConfs



class Optimizer:

    def __init__(self):
        self.workdir_name = "supra_running_opt"
        self.coord_file_name = "coord"
        self.temp_file_name = "temp.xyz"
        self.work_struc_name = "struc.xyz"
        self.opt_struc_name = "opt_struc.xyz"
    


    def optimize_structure_uff(self, coords: dict, n: int = None):
        xyz_string = f"{len(coords.keys())}\n\n"
        for atom, (x, y, z) in coords.items():
            xyz_string = xyz_string + f"{get_element(atom)}\t{x}\t{y}\t{z}\n"

        mol = Chem.MolFromXYZBlock(xyz_string)
        mol = Chem.Mol(mol)
        rdDetermineBonds.DetermineBonds(mol)

        #res = UFFOptimizeMoleculeConfs(mol, maxIters=1000)
        res = MMFFOptimizeMoleculeConfs(mol, maxIters=1000)
        energy = res[0][1]

        opt_struc_file = os.path.join(os.getcwd(), f"conformer{n}.xyz")
        with open(opt_struc_file, "w") as outfile:
            print(len(coords.keys()), file=outfile)
            print(f"UFF-Energy = {energy}", file=outfile)
            for i, atom in enumerate(mol.GetAtoms()):
                pos = mol.GetConformer().GetAtomPosition(i)
                print(f"{atom.GetSymbol()}\t{pos.x}\t{pos.y}\t{pos.z}", file=outfile)



    """
    def optimize_structure_uff(self, coords: dict, n: int=None):
        workdir = os.path.join(os.getcwd(), f"{self.workdir_name}{n}")
        workdir = os.path.abspath(workdir)
        os.makedirs(workdir)

        new_xyz_file = os.path.join(workdir, "struc.xyz")
        os.system(f"touch {new_xyz_file}")
        with open(new_xyz_file, "w") as struc_to_optimize:
            number_of_atoms = len(coords.keys())
            print(number_of_atoms, file=struc_to_optimize, end="\n\n")
            for atom, (x, y, z) in coords.items():
                print(f"{get_element(atom)}\t{x}\t{y}\t{z}", file=struc_to_optimize)
        
        control_file = os.path.join(workdir, "control")
        coord_file = os.path.join(workdir, "coord")
        os.system(f"touch {control_file}")
        with open(control_file, "w") as control:
            print("$symmetry c1", file=control)
            print("$uff", file=control)
            print("      2500         1          0 ! maxcycle,modus,nqeq", file=control)
            print("    111111                      ! iterm", file=control)
            print("  0.10D-07  0.10D-04            ! econv,gconv", file=control)
            print("      0.00  1.10                ! qtot,dfac", file=control)
            print("  0.10D+03  0.10D-04       0.30 ! epssteep,epssearch,dqmax", file=control)
            print("        25      0.10       0.00 ! mxls,dhls,ahls", file=control)
            print("      1.00      0.00       0.00 ! alpha,beta,gamma", file=control)
            print("         F         F          F ! transform,lnumhess,lmd", file=control)
            print("$end", file=control)
        
        with open(coord_file, "w") as f:
            subprocess.run(
                args=["x2t", self.work_struc_name, ">", self.coord_file_name],
                cwd=workdir,
                stdout=f,
                stderr=subprocess.DEVNULL
            )
        
        subprocess.run(
            args=["uff"],
            cwd=workdir,
            stdout=subprocess.DEVNULL,
            stderr=subprocess.DEVNULL
        )

        with open(os.path.join(workdir, "uffenergy"), "r") as infile:
            energy = float(infile.read().split()[-2])

        opt_struc = os.path.join(workdir, self.opt_struc_name)
        temp_file = os.path.join(workdir, self.temp_file_name)
        with open(temp_file, "w") as f:
            subprocess.run(
                args=["t2x", self.coord_file_name, ">", self.temp_file_name],
                cwd=workdir,
                stdout=f,
                stderr=subprocess.DEVNULL
            )
        
        with open(temp_file, "r") as infile:
            with open(opt_struc, "w") as outfile:
                for i, line in enumerate(infile):
                    if i == 1:
                        print(f"Energy = {energy}", file=outfile)
                    else:
                        print(line, file=outfile, end="")
        
        os.system(
            f"mv {opt_struc} conformer{n}.xyz ; \
              rm -r {workdir}"
        )
    """

    

    def optimize_structure_xtb(self, struc_folder: str, struc_file: str, chrg: int=None, n: int=None):
        workdir = os.path.join(struc_folder, self.workdir_name+str(n))
        xtb_args = ["xtb", "--opt", "normal"]
        if chrg:
            xtb_args.append("--chrg")
            xtb_args.append(str(chrg))
        os.system(
            f"mkdir {workdir} ; \
              mv {os.path.join(struc_folder, struc_file)} {workdir}"
        )
        working_struc = os.path.join(workdir, struc_file)
        xtb_args.append(working_struc)
        subprocess.run(
            args=xtb_args, 
            cwd=workdir, 
            stdout=subprocess.DEVNULL, 
            stderr=subprocess.DEVNULL
        )
        xtb_struc_file = os.path.join(workdir, "xtbopt.xyz")
        opt_struc_file = os.path.join(struc_folder, f"conformer{n}.xyz")
        os.system(f"mv {xtb_struc_file} {opt_struc_file}")
        os.system(f"rm -rf {workdir}")



    def refine_structures_xtb(self, path_to_strucs: str, chrg: int=None):
        strucs_list = os.listdir(path_to_strucs)
        print("Performing refining optimizations...")
        for i, struc in enumerate(strucs_list):
            self.optimize_structure_xtb(
                struc_folder=path_to_strucs,
                struc_file=struc,
                chrg=chrg, 
                n=i
            )
        print("Refining optimizations done.")
