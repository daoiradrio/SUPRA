import os
import subprocess



class Optimizer:

    def __init__(self):
        self.workdir_name = "supra_running_opt"

    
    def optimize_structure_uff():
        pass

    
    def optimize_structure_xtb(struc_file: str, chrg: int=None, n: int=None):
        xtb_args = ["xtb", "--opt", "normal"]
        if chrg:
            xtb_args.append("--chrg")
            xtb.args.append(str(chrg))
        os.system(
            f"mkdir {self.workdir_name} ; mv {struc_file} {self.workdir_name}"
        )
        struc_path = os.path.join(
            os.path.abspath(self.workdir_name),
            struc_file
        )
        subprocess.run(
            args=["xtb", "--opt", "normal"]+[struc_path], 
            cwd=new_workdir, 
            stdout=subprocess.DEVNULL, 
            stderr=subprocess.DEVNULL
        )
        xtb_struc_file = os.path.join(self.workdir_name, "xtbopt.xyz")
        opt_struc_file = f"../conformer{n}.xyz"
        os.system(f"mv {xtb_struc_file} {opt_struc_file}")
        os.system(f"rm -rf {self.workdir_name}")


    def refine_structures_xtb():
        strucs_list = os.listdir(path_to_strucs)
        print("Performing refining optimizations...")
        for i, struc in enumerate(strucs_list):
            self.optimize_structure_xtb(struc, chrg, i)
        print("Refining optimizations done.")
