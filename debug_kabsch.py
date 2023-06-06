import os
from utils.analyzer import Analyzer
from SUPRAConformer.structure import Structure

file1 = os.sys.argv[1]
file2 = os.sys.argv[2]
analyzer = Analyzer()
struc1 = Structure(file1)
struc2 = Structure(file2)
analyzer.doubles(struc1.coords, struc2.coords)
