import numpy as np


E = 200e9
L = 10.0
D = 0.1
R = D / 2
I = np.pi / 4 * R**4
M = 2 * np.pi * E * I / L
print("M\n", M)
