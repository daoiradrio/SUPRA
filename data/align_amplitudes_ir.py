file1 = "/home/baum/SUPRA/data/Dehydroquinidin_ignore_methyl_120/ir/normalized_scaled_ir.dat"
file2 = "/home/baum/SUPRA/data/Dehydroquinidin_ignore_methyl_120/ir/normalized_sim_ir.dat"

new_file1 = "/home/baum/SUPRA/data/Dehydroquinidin_ignore_methyl_120/ir/scaled_aligned_ir.dat"
new_file2 = "/home/baum/SUPRA/data/Dehydroquinidin_ignore_methyl_120/ir/aligned_sim_ir.dat"

min1 = 0.00000000
min2 = 0.07153081
xs1 = []
ys1 = []
xs2 = []
ys2 = []

min1 = 1 - min1
min2 = 1 - min2

with open(file1, "r") as infile:
    for i, line in enumerate(infile):
        data = line.split()
        xs1.append(float(data[0]))
        ys1.append(float(data[1]))

ys1 = [1-y for y in ys1]

with open(new_file1, "w") as outfile:
    for x, y in zip(xs1, ys1):
        y = 1 - y / min1
        print(f"{x:15.8f}\t{y:15.8f}", file=outfile)

with open(file2, "r") as infile:
    for i, line in enumerate(infile):
        data = line.split()
        xs2.append(float(data[0]))
        ys2.append(float(data[1]))

ys2 = [1-y for y in ys2]

with open(new_file2, "w") as outfile:
    for x, y in zip(xs2, ys2):
        y = 1 - y / min2
        print(f"{x:15.8f}\t{y:15.8f}", file=outfile)
