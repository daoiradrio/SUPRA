import os

from utils.symmetry import Symmetry
from utils.rotationaxis import RotationAxis
from SUPRAConformer.structure import Structure



files = "testcases/"
try:
    f = open(os.path.join(files, "Alanin.xyz"))
    f.close()
except:
    files = "tests/testcases/"