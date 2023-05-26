import os
import subprocess



class Optimizer:

    def __init__(self):
        self.workdir_name = "supra_running_opt"

    
    def optimize_structure_uff():
        pass

    
    def optimize_structure_xtb(self, struc_folder: str, struc_file: str, chrg: int=None, n: int=None):
        workdir = os.path.join(struc_folder, self.workdir_name+str(n))
        xtb_args = ["xtb", "--opt", "normal"]
        if chrg:
            xtb_args.append("--chrg")
            xtb.args.append(str(chrg))
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
