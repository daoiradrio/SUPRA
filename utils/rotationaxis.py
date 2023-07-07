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