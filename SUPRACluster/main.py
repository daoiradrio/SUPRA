from clusterstructure import ClusterStructure
from clustergenerator import ClusterGenerator
import numpy as np



filename = "../inputfiles/Alanin.xyz"
struc = ClusterStructure()
struc.get_structure(filename)
struc.find_hbs()
gen = ClusterGenerator(struc)
gen.add_monomer_at_don(struc, "N8", "H14")
