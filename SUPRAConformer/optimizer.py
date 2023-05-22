import os
import subprocess



class Optimizer:

    def __init__(self):
        self.workdir_name = "supra_running_opt"


    def refine_ff_opts(self, path_to_strucs: str, chrg: int=None) -> None:
        strucs_list = os.listdir(path_to_strucs)
        workdir = os.path.join(os.path.abspath(path_to_strucs), self.workdir_name)
        xtb_args = ["xtb", "--opt", "normal"]
        if chrg:
            xtb_args.append("--chrg")
            xtb_args.append(str(chrg))
        print("Performing refining optimizations...")
        for i, struc in enumerate(strucs_list):
            new_workdir = f"{workdir}{i}"
            os.system(f"mkdir {new_workdir} ; mv {os.path.join(path_to_strucs, struc)} {new_workdir}")
            struc_path = os.path.join(os.path.abspath(new_workdir), struc)
            subprocess.run(args=xtb_args+[struc_path], cwd=new_workdir, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)
            xtb_struc_file = os.path.join(new_workdir, "xtbopt.xyz")
            opt_struc_file = os.path.join(path_to_strucs, f"conformer{i}.xyz")
            os.system(f"mv {xtb_struc_file} {opt_struc_file}")
            os.system(f"rm -rf {new_workdir}")
