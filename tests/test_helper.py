import os

from utils.helper import get_element, get_number




files = "testcases/"
try:
    f = open(os.path.join(files, "Alanin.xyz"))
    f.close()
except:
    files = "tests/testcases/"



def test_get_element():
    atom1 = "C0"
    atom2 = "N11"
    atom3 = "Cl7"
    atom4 = "Br42"

    check1 = (get_element(atom1) == "C")
    check2 = (get_element(atom2) == "N")
    check3 = (get_element(atom3) == "Cl")
    check4 = (get_element(atom4) == "Br")

    assert (check1 and check2 and check3 and check4)



def test_get_number():
    atom1 = "C0"
    atom2 = "N11"
    atom3 = "Cl7"
    atom4 = "Br42"

    check1 = (get_number(atom1) == 0)
    check2 = (get_number(atom2) == 11)
    check3 = (get_number(atom3) == 7)
    check4 = (get_number(atom4) == 42)

    assert (check1 and check2 and check3 and check4)
