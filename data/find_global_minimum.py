import os
from SUPRAConformer.structure import Structure



ensemble_path = "/home/baum/SUPRA/data/Butan-2-ol_all_groups_60/ensemble"
conformer_files = os.listdir(ensemble_path)

E_min = 1000000.0
conf_min = None
for file in conformer_files:
    struc = Structure()
    struc.read_xyz(
        filename=os.path.join(ensemble_path, file),
        read_energy=True
    )
    if struc.energy < E_min:
        E_min = struc.energy
        conf_min = file.split(".")[0]

print(conf_min)