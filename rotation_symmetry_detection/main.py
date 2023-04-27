import os
import numpy as np
from SUPRAConformer.structure import Structure



tolerance_primary = 5e-2



def geometric_center(coords: list):
    center = np.array([0.0, 0.0, 0.0])
    for coord in coords:
        center += coord
    n = len(coords)
    center = center / n
    return center 


def find_c2(atoms, coords, geometric_center, distances_from_center):
    global tolerance_primary
    center = np.array([0.0, 0.0, 0.0])
    r = 0.0
    # loop over all atoms
    for i, atom1 in enumerate(atoms):
        # loop over all atoms up to the previous atom
        for j, atom2 in enumerate(atoms[:i]):
            # continue with next iteration if elements don't match
            if Structure.get_element(atom1) != Structure.get_element(atom2):
                continue
            # if both atoms are not equally close to the geometric center (within some range of tolerance)
            # then ???, continue with next iteration
            if abs(distances_from_center[i] - distances_from_center[j] > tolerance_primary):
                continue
            # center of vector between the two current atoms
            center = 0.5 * (coords[i] - coords[j])
            # distance of center of vector to geometric center
            r = (center - geometric_center)**2
            r = np.sqrt(r)
            if r < 5 * tolerance_primary:
                pass




path = "/home/baum/SUPRA/inputfiles/"
xyz_file = "Tyrosin.xyz"
xyz_path = os.path.join(path, xyz_file)
mol = Structure(xyz_path)
atoms = [atom for atom in mol.coords.keys()]
coords = [coord for coord in mol.coords.values()]
center = geometric_center(coords)
distances_from_center = [np.linalg.norm(center - coord) for coord in coords]
find_c2(atoms, coords, center, distances_from_center)
