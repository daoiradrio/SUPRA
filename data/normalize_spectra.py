import os



filepath = "/home/dario/SUPRA/data/Serin_all_groups_120/weighted_ir.dat"

xs = []
ys = []
with open(filepath, "r") as infile:
    for i, line in enumerate(infile):
        data = line.split()
        xs.append(float(data[0]))
        ys.append(float(data[1]))

max_y = max(ys)
for i, y in enumerate(ys):
    ys[i] = 1 - y / max_y

with open("normalized_inverted_ir.dat", "w") as outfile:
    for x, y in zip(xs, ys):
        print(f"{x:15.8f}\t{y:15.8f}", file=outfile)
