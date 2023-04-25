from clusterstructure import ClusterStructure



def output(structure: ClusterStructure, new_coords: list):
    with open("out.xyz", "w") as outfile:
        n = len(structure.coords.keys()) + len(new_coords)
        print(n, file=outfile, end="\n\n")
        for atom, coords in structure.coords.items():
            element = structure.get_element(atom)
            x = coords[0]
            y = coords[1]
            z = coords[2]
            print(f"{element}\t{x}\t{y}\t{z}", file=outfile)
        for coords in new_coords:
            x = coords[0]
            y = coords[1]
            z = coords[2]
            print(f"F\t{x}\t{y}\t{z}", file=outfile)
    


filename = "../inputfiles/Alanin.xyz"
struc = ClusterStructure()
struc.get_structure(filename)
struc.find_hbs()
new_coords = [acc_vec for acc_vec in struc.hb_acc_vec.values()]
output(struc, new_coords)
