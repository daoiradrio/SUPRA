from clusterstructure import ClusterStructure
from clustergenerator import ClusterGenerator



filename = "../inputfiles/Alanin.xyz"
struc = ClusterStructure()
struc.get_structure(filename)
struc.find_hbs()
gen = ClusterGenerator(struc)
#gen.add_monomer(struc, "N8", "H14")
