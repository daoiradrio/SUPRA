import os

import numpy as np

from utils.rotationaxis import RotationAxis



files = "testcases/"
try:
    f = open(os.path.join(files, "Alanin.xyz"))
    f.close()
except:
    files = "tests/testcases/"



def test_rotate_atom():
    tol = 0.00001

    from_coords = np.array([0.0, 0.0, 0.0])
    to_coords = np.array([0.0, 0.0, 1.0])
    coords = np.array([1.0, 0.0, 0.0])
    deg = 90.0

    new_coords1 = RotationAxis.rotate_atom(
        from_coords=from_coords,
        to_coords=to_coords,
        coords=coords,
        deg=deg
    )

    from_coords = np.array([0.0, 0.0, 0.0])
    to_coords = np.array([1.0, 1.0, 1.0])
    coords = np.array([1.0, 1.0, 0.0])
    deg = 180.0

    new_coords2 = RotationAxis.rotate_atom(
        from_coords=from_coords,
        to_coords=to_coords,
        coords=coords,
        deg=deg
    )

    assert (
        abs(new_coords1[0] - 0.0) < tol and
        abs(new_coords1[1] - 1.0) < tol and
        abs(new_coords1[2] - 0.0) < tol and
        abs(new_coords2[0] - 1/3) < tol and
        abs(new_coords2[1] - 1/3) < tol and
        abs(new_coords2[2] - 4/3) < tol
    )
