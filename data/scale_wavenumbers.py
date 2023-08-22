filepath = "/home/baum/SUPRA/data/Dehydroquinidin_ignore_methyl_120/ir/normalized_glob_min_ir.dat"
new_file = "/home/baum/SUPRA/data/Dehydroquinidin_ignore_methyl_120/ir/normalized_scaled_glob_min_ir.dat"

xs = []
ys = []
with open(filepath, "r") as infile:
    for i, line in enumerate(infile):
        data = line.split()
        xs.append(float(data[0]))
        ys.append(float(data[1]))
with open(new_file, "w") as outfile:
    for x, y in zip(xs, ys):
        x = 0.97 * x
        print(f"{x:15.8f}\t{y:15.8f}", file=outfile)
