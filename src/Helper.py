import numpy as np
import os
import shutil



# covalence radii of chain atoms for single bonds, double bonds and triple bonds in angstrom
# for calculating the connectivity
#
# Sources:
# single bonds: P. Pyykkö, M. Atsumi,Chemistry – A European Journal2009,15, 186–197
# double bonds: P. Pyykkö, M. Atsumi,Chemistry – A European Journal2009,15, 12770–12779
# triple bonds: P. Pyykkö, S. Riedel, M. Patzschke,Chemistry – A European Journal2005,11, 3511–3520
covalence_radii = {"C": [0.75, 0.67, 0.60], "N": [0.71, 0.60, 0.54], "O": [0.66, 0.57, 0.53], "H": [0.32],
                   "B": [0.85, 0.78, 0.73], "F": [0.64], "Cl": [0.99], "Br": [1.14], "I": [1.33]}


# maximum valences of chain atoms for regulation of loop cycles in building self.structure (Structure-class)
valences = {"C": 4, "N": 3, "O": 2, "H": 1, "B": 3, "F": 1, "Cl": 1, "Br": 1, "I": 1}


#WBB-LÄNGEN FÜR VERSCHIEDENE ELEMENTE EINPFLEGEN (SERIÖSE QUELLE SUCHEN)
# hb_lengths = {}



# translating of an atom labels in element symbols
def get_element(label: str) -> str:
    # if letter on second index: two letter element symbol, return first two chars of label
    if label[1].isalpha():
        return (label[0] + label[1])
    # if no letter on second index: one letter element symbol, return first char of label
    else:
        return label[0]


def get_number(label: str) -> int:
    if label[-2].isdigit():
        return int(label[-2] + label[-1])
    else:
        return int(label[-1])


#--> ABFRAGE HINZUFÜGEN OB WINKEL BEREITS IN BOGENMAß, FALLS NICHT SELBST UMWANDELN
#    WAHRSCHEINLICH ÜBER WERTEBEREICH MACHBAR <-- NÖ, ÜBERGABE PARAMETER rad als bool
def rotation(vec2, vec1, angle: float, coords, rad: bool = False) -> np.array:
    if type(vec2) == list:
        vec2 = np.array(vec2)
    if type(vec1) == list:
        vec1 = np.array(vec1)
    if type(coords) == list:
        coords = np.array(coords)
    if not rad:
        angle = np.deg2rad(angle)

    axis = vec2 - vec1
    axis = axis / np.sqrt(np.dot(axis, axis))
    coords = coords - vec1

    coords = np.dot(axis, np.dot(axis, coords)) \
             + np.cos(angle) * np.cross(np.cross(axis, coords), axis) \
             + np.sin(angle) * np.cross(axis, coords)

    coords = coords + vec1

    return coords


def is_hb_don(label: str) -> bool:
    flag = False
    element = get_element(label)
    if element == "O" or element == "N":
        flag = True
    return flag


def is_hb_acc(label: str) -> bool:
    flag = False
    if get_element(label) == "H":
        flag = True
    return flag


def generate_folder(foldername: str) -> str:
    abs_path_new_folder = os.path.abspath(foldername)
    if os.getcwb != abs_path_new_folder:
        if os.path.exists(abs_path_new_folder):
            shutil.rmtree(abs_path_new_folder, ignore_errors=True)
        os.makedirs(abs_path_new_folder)
    return abs_path_new_folder


def guess_molecule_name(infile: str, molecule_name: str) -> str:
    if molecule_name == None:
        molecule_name = os.path.splittext(os.path.basename(infile))[0]
    return molecule_name
