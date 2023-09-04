file1 = "/home/baum/SUPRA/data/Dehydroquinidin_ignore_methyl_120/new_vcd/normalized_vcd.dat"
file2 = "/home/baum/SUPRA/data/Dehydroquinidin_ignore_methyl_120/new_vcd/sim_vcd.dat"

new_file1 = "/home/baum/SUPRA/data/Dehydroquinidin_ignore_methyl_120/new_vcd/aligned_vcd.dat"
new_file2 = "/home/baum/SUPRA/data/Dehydroquinidin_ignore_methyl_120/new_vcd/final_sim_vcd.dat"

max1 = 1.00000000
max2 = 0.8230946684894054
xs1 = []
ys1 = []
xs2 = []
ys2 = []

with open(file1, "r") as infile:
    for i, line in enumerate(infile):
        data = line.split()
        xs1.append(float(data[0]))
        ys1.append(float(data[1]))

with open(new_file1, "w") as outfile:
    for x, y in zip(xs1, ys1):
        y = y / abs(max1)
        print(f"{x:15.8f}\t{y:15.8f}", file=outfile)

with open(file2, "r") as infile:
    for i, line in enumerate(infile):
        data = line.split()
        xs2.append(float(data[0]))
        ys2.append(float(data[1]))

with open(new_file2, "w") as outfile:
    for x, y in zip(xs2, ys2):
        y = y / abs(max2)
        print(f"{x:15.8f}\t{y:15.8f}", file=outfile)
