import numpy as np

from scipy.spatial.transform import Rotation



axis1 = np.array([1.0, 0.5, 1.0])
axis2 = np.array([-0.5, 1.0, -1.0])
axis = (axis2 - axis1) / np.linalg.norm(axis2 - axis1)
p = [1, 0, 0]
R = Rotation.from_rotvec(np.pi/2 * axis)

print(R.apply(p))
print(p)