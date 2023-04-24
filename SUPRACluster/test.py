import sys
import numpy as np



angle = float(sys.argv[1])
arccos = np.arccos(angle)
print(np.rad2deg(arccos))
