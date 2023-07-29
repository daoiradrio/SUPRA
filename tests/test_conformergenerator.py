import os

from SUPRAConformer.structure import Structure
from SUPRAConformer.conformergenerator import ConformerGenerator



files = "testcases/"
try:
    f = open(os.path.join(files, "Alanin.xyz"))
    f.close()
except:
    files = "tests/testcases/"



def test_find_torsions():
    assert True

def test_find_cylces():
    assert True

def test_peptidebonds():
    assert True

def test_clashes():
    assert True