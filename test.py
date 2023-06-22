with open("Alanin.xyz", "r") as f:
    for i, line in enumerate(f):
        if i == 1:
            print(line.split()[-1])