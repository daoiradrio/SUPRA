filepath = "/home/dario/SUPRA/data/Dehydroquinidin_ignore_methyl_120/weighted_ir.dat"

xs = []
ys = []
with open(filepath, "r") as infile:
    for i, line in enumerate(infile):
        data = line.split()
        xs.append(float(data[0]))
        ys.append(float(data[1]))

max_y = abs(max(ys))
min_y = abs(min(ys))
if (max_y > min_y):
    ref = max_y
else:
    ref = min_y
for i, y in enumerate(ys):
    ys[i] = 1 - y / ref
    #ys[i] = y / ref

with open("normalized_ir.dat", "w") as outfile:
    for x, y in zip(xs, ys):
        print(f"{x:15.8f}\t{y:15.8f}", file=outfile)
