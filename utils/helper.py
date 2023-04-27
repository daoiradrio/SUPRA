import numpy as np



# covalence radii of chain atoms for single bonds, double bonds and triple bonds in angstrom
# for calculating the connectivity
#
# Sources:
# single bonds: P. Pyykkö, M. Atsumi,Chemistry – A European Journal2009,15, 186–197
# double bonds: P. Pyykkö, M. Atsumi,Chemistry – A European Journal2009,15, 12770–12779
# triple bonds: P. Pyykkö, S. Riedel, M. Patzschke,Chemistry – A European Journal2005,11, 3511–3520
covalence_radii = {
    "C": [0.75, 0.67, 0.60], "N": [0.71, 0.60, 0.54], "O": [0.66, 0.57, 0.53], "H": [0.32], "B": [0.85, 0.78, 0.73],
    "F": [0.64], "Cl": [0.99], "Br": [1.14], "I": [1.33]
}
covalence_radii_single = {
    "C": 0.75, "N": 0.71, "O": 0.66, "H": 0.32, "B": 0.85, "F": 0.64, "Cl": 0.99, "Br": 1.14, "I": 1.33
}
covalence_radii_double = {
    "C": 0.67, "N": 0.60, "O": 0.57, "B": 0.78
}
covalence_radii_triple = {
    "C": 0.60, "N": 0.54, "O": 0.53, "B": 0.73
}


# maximum valences of chain atoms for regulation of loop cycles in building self.structure (Structure-class)
valences = {
    "C": 4, "N": 3, "O": 2, "H": 1, "B": 3, "F": 1, "Cl": 1, "Br": 1, "I": 1
}


atomic_masses_symbols = {
    "H": 1.00797, "He": 4.00260,
    "Li": 6.941, "Be": 9.01218, "B": 10.81, "C": 12.011, "N": 14.0067, "O": 15.9994, "F": 18.998403, "Ne": 0.179,
}

atomic_masses_numbers = {
    1: 1.00797, 2: 4.00260,
    3: 6.941, 4: 9.01218, 5: 10.81, 6: 12.011, 7: 14.0067, 8: 15.9994, 9: 18.998403, 10: 0.179,
}


# combinations of angle increments
increment_combinations = {
    30: [30, 45], 45: [45, 60], 60: [60, 90], 90: [90, 120], 120: [120, 180], 180: [180]
}


#WBB-LÄNGEN FÜR VERSCHIEDENE ELEMENTE EINPFLEGEN (SERIÖSE QUELLE SUCHEN)
# hb_lengths = {}



def rotation(coord_rotation_atom: np.array, coord_axis_atom1: np.array, coord_axis_atom2: np.array, rotation_angle: float, deg: bool=False) -> np.array:
    if deg:
        rotation_angle = np.deg2rad(rotation_angle)
    axis = coord_axis_atom2 - coord_axis_atom1
    axis = axis / np.linalg.norm(axis)
    coord_rotation_atom = coord_rotation_atom - coord_axis_atom1
    coord_rotation_atom = np.dot(axis, np.dot(axis, coord_rotation_atom)) \
             + np.cos(rotation_angle) * np.cross(np.cross(axis, coord_rotation_atom), axis) \
             + np.sin(rotation_angle) * np.cross(axis, coord_rotation_atom)
    coord_rotation_atom = coord_rotation_atom + coord_axis_atom1
    return coord_rotation_atom


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


def is_hb_don(label: str) -> bool:
    flag = False
    element = get_element(label)
    if element in ["O", "N"]: # WELCHE ELEMENTE NOCH?? HALOGENE??
        flag = True
    return flag


def is_hb_acc(label: str) -> bool:
    flag = False
    if get_element(label) == "H":
        flag = True
    return flag


def geometric_center(coords: dict) -> np.array:
    center = np.array([0.0, 0.0, 0.0])
    for coord in coords.values():
        center += coord
    n = len(coords)
    center = center / n
    return center


def mass_center(coords: dict) -> np.array:
    sum_m = 0.0
    com = np.array([0.0, 0.0, 0.0])
    for atom, coord in coords.items():
        element = get_element(atom)
        m = atomic_masses_symbols[element]
        com += m * coord
        sum_m += m
    com *= (1/m)
    return com


def inertia_tensor(coords: dict) -> np.array:
    com = mass_center(coords)
    I = np.zeros((3,3))
    for i in range(3):
        for j in range(3):
            for atom, coord in coords.items():
                r = coord - com
                element = get_element(atom)
                m = atomic_masses_symbols[element]
                I[i][j] -= (m * r[i]*r[j])
                if i == j:
                    I[i][j] += (m * np.dot(r, r))
    return I