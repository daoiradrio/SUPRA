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
    def new_rotate_atom(p, x1, x2, theta):
        p = np.append(p, 1)
        p = [[pp] for pp in p]
        x1, y1, z1 = x1
        x2, y2, z2 = x2

        U = [x2-x1, y2-y1, z2-z1]
        U = np.array(U) / np.sqrt(np.dot(U,U))
        a,b,c = U
        d = np.sqrt(b**2 + c**2)

        T = [
            [1,0,0,-x1],
            [0,1,0,-y1],
            [0,0,1,-z1],
            [0,0,0,1  ]
        ]
        T_inv = [
            [1,0,0,x1],
            [0,1,0,y1],
            [0,0,1,z1],
            [0,0,0,1 ]
        ]

        R_x = [
            [1,0,0,0     ],
            [0,c/d,-b/d,0],
            [0,b/d,c/d,0 ],
            [0,0,0,1     ]
        ]
        R_x_inv = [
            [1,0,0,0     ],
            [0,c/d,b/d,0 ],
            [0,-b/d,c/d,0],
            [0,0,0,1     ]
        ]

        R_y = [
            [d,0,-a,0],
            [0,1,0,0 ],
            [a,0,d,0 ],
            [0,0,0,1 ]
        ]
        R_y_inv = [
            [d,0,a,0 ],
            [0,1,0,0 ],
            [-a,0,d,0],
            [0,0,0,1 ]
        ]

        ct = np.cos(theta)
        st = np.sin(theta)
        R_z = [
            [ct,st,0,0 ],
            [-st,ct,0,0],
            [0,0,1,0   ],
            [0,0,0,1   ]
        ]

        p2 = np.dot(T, p)
        p2 = np.dot(R_x, p2)
        p2 = np.dot(R_y, p2)
        p2 = np.dot(R_z, p2)
        p2 = np.dot(R_y_inv, p2)
        p2 = np.dot(R_x_inv, p2)
        p2 = np.dot(T_inv, p)

        return p2[0][:3]



    @staticmethod
    def matrix_multiply(*matrices):
        if len(matrices) == 1:
            return matrices
        else:
            try:
                m_other = RotationAxis.matrix_multiply(*matrices[1:])
                return np.matmul(matrices[0], m_other)
            except:
                #print(matrices[0])
                #print(m_other)
                raise