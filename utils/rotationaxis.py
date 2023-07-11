import numpy as np



class RotationAxis:

    def __init__(self, from_coords: np.array, to_coords: np.array):
        self.from_coords = from_coords
        self.to_coords   = to_coords
        self.axis        = (to_coords - from_coords) / np.linalg.norm(to_coords - from_coords)


    

    def rotate(self, coords: np.array, deg: float) -> np.array:
        angle = np.deg2rad(deg)

        new_coords = coords - self.from_coords
        new_coords = np.dot(self.axis, np.dot(self.axis, new_coords)) \
                     + np.cos(angle) * np.cross(np.cross(self.axis, new_coords), self.axis) \
                     + np.sin(angle) * np.cross(self.axis, new_coords)
        new_coords = new_coords + self.from_coords

        return new_coords



    @staticmethod
    def rotate_atom(from_coords: np.array, to_coords: np.array, coords: np.array, deg: float) -> np.array:
        axis = (to_coords - from_coords) / np.linalg.norm(to_coords - from_coords)
        angle = np.deg2rad(deg)

        new_coords = coords - from_coords
        new_coords = np.dot(axis, np.dot(axis, new_coords)) \
                     + np.cos(angle) * np.cross(np.cross(axis, new_coords), axis) \
                     + np.sin(angle) * np.cross(axis, new_coords)
        new_coords = new_coords + from_coords

        return new_coords



    @staticmethod
    def new_rotate_atom(p, x1, x2, angle):
        a = (x2 - x1) / np.linalg.norm(x2 - x1)
        ax, ay, az = a
        theta = np.deg2rad(angle)

        R = np.zeros((3, 3))
        R[0][0] = ax**2 + (ay**2 + az**2) * np.cos(theta)
        R[0][1] = ax * ay * (1 - np.cos(theta)) - az * np.sin(theta)
        R[0][2] = ax * az * (1 - np.cos(theta)) + ay * np.sin(theta)
        R[1][0] = ax * ay * (1 - np.cos(theta)) + az * np.sin(theta)
        R[1][1] = ay**2 + (ax**2 + az**2) * np.cos(theta)
        R[1][2] = ay * az * (1 - np.cos(theta)) - ax * np.sin(theta)
        R[2][0] = ax * az * (1 - np.cos(theta)) - ay * np.sin(theta)
        R[2][1] = ay * az * (1 - np.cos(theta)) + ax * np.sin(theta)
        R[2][2] = az**2 + (ax**2 + ay**2) * np.cos(theta)

        new_p = p - x1
        new_p = np.dot(R, new_p)
        new_p = new_p + x1

        return new_p