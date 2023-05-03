from clusterstructure import ClusterStructure
from clustergenerator import ClusterGenerator
from copy import deepcopy



def output(structure: ClusterStructure):
    with open("out.xyz", "w") as outfile:
        n = len(structure.coords.keys())
        print(n, file=outfile, end="\n\n")
        for atom, coords in structure.coords.items():
            element = structure.get_element(atom)
            x = coords[0]
            y = coords[1]
            z = coords[2]
            print(f"{element}\t{x}\t{y}\t{z}", file=outfile)
    


filename = "../inputfiles/Alanin.xyz"
struc = ClusterStructure()
struc.get_structure(filename)
gen = ClusterGenerator(struc)
copy_struc = deepcopy(struc)
gen.add_monomer(struc, "H14", copy_struc, "N8")
output(struc)