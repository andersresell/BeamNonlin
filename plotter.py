import numpy as np
import matplotlib.pyplot as plt
import os

# fmt: off

if __name__ == "__main__":
    os.chdir("testing/output")

    data = np.genfromtxt("0.csv", delimiter=",", skip_header=1)
    N = data[0, 0]
    n_steps = data[0, 3]
    n_write = data[0, 4]

    for n in range(n_write):
        if n % n_write != 0: continue

        data = np.genfromtxt(str(n) + ".csv", delimiter=",", skip_header=4)
        print(data)
        break
