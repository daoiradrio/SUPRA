import os

import numpy as np

from scipy import constants
from SUPRAConformer.structure import Structure

Eh_to_J = 4.3597447222071e-18



def boltzmann(E: float):
    return np.exp(-Eh_to_J*E/(273.15 * constants.Boltzmann))

def get_boltzmann_factors(energies: list):
    E_min = min(energies)
    normalized_energies = [E-E_min for E in energies]
    tmp = []
    for E in normalized_energies:
        tmp.append(boltzmann(E))
    Z = sum(tmp)
    return [(1/Z)*val for val in tmp]



path = "/home/dario/SUPRA/tmp"
molset = "Tyrosin_all_groups_90"
#folder = os.path.join(path, molset)
folder = path

conformer_files = os.listdir(os.path.join(folder, "ensemble"))
ir_files = os.listdir(os.path.join(folder, "ir"))
vcd_files = os.listdir(os.path.join(folder, "vcd"))

energies = []
for conformer_file in conformer_files:
    new_struc = Structure()
    new_struc.read_xyz(filename=os.path.join(folder, "ensemble", conformer_file), read_energy=True)
    energies.append(new_struc.energy)
boltzmann_factors = get_boltzmann_factors(energies)

vcd_xs = []
with open(os.path.join(folder, "vcd", vcd_files[0]), "r") as infile:
    for i, line in enumerate(infile):
        data = line.split()
        if i == 0:
            continue
        vcd_xs.append(float(data[0]))

weighted_vcd_ys = [0.0 for _ in vcd_xs]
for boltzmann_factor, vcd_file in zip(boltzmann_factors, vcd_files):
    with open(os.path.join(folder, "vcd", vcd_file), "r") as infile:
        for i, line in enumerate(infile, start=-1):
            data = line.split()
            if i == -1:
                continue
            weighted_vcd_ys[i] += boltzmann_factor*float(data[1])

ir_xs = []
with open(os.path.join(folder, "ir", ir_files[0]), "r") as infile:
    for i, line in enumerate(infile):
        data = line.split()
        if i == 0:
            continue
        ir_xs.append(float(data[0]))

weighted_ir_ys = [0.0 for _ in ir_xs]
for boltzmann_factor, ir_file in zip(boltzmann_factors, ir_files):
    with open(os.path.join(folder, "ir", ir_file), "r") as infile:
        for i, line in enumerate(infile, start=-1):
            data = line.split()
            if i == -1:
                continue
            weighted_ir_ys[i] += boltzmann_factor*float(data[1])

new_file = "weighted_ir.dat"
with open(new_file, "w") as outfile:
    for x, y in zip(ir_xs, weighted_ir_ys):
        print(f"{x:15.8f}\t{y:15.8f}", file=outfile)

new_file = "weighted_vcd.dat"
with open(new_file, "w") as outfile:
    for x, y in zip(vcd_xs, weighted_vcd_ys):
        print(f"{x:15.8f}\t{y:15.8f}", file=outfile)