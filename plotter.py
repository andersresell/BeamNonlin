import numpy as np
import matplotlib.pyplot as plt
import os
from mpl_toolkits.mplot3d import Axes3D

# fmt: off

#need this to make ctrl+c close plot window
import signal
signal.signal(signal.SIGINT, signal.SIG_DFL)

SMALL_VAL = 1e-6

def plot_triad(ax, x,U):
    assert x.shape==(3,)
    assert U.shape==(3,3)
    ax.quiver(x[0],x[1],x[2],U[0,0],U[0,1],U[0,2,],color="r")
    ax.quiver(x[0],x[1],x[2],U[1,0],U[1,1],U[1,2,],color="g")
    ax.quiver(x[0],x[1],x[2],U[2,0],U[2,1],U[2,2,],color="b")


if __name__ == "__main__":
    os.chdir("testing/output")

    header = np.genfromtxt("0.csv", delimiter=",", skip_header=1,max_rows=1)
    N = int(header[0])
    n_steps = int(header[1])
    n_write = int(header[2])
       
    ax = plt.figure().add_subplot(111,projection="3d")
    for n in range(n_steps):
        if n % n_write != 0: continue

        header = np.genfromtxt(str(n)+".csv", delimiter=",", skip_header=1,max_rows=1)
        t = header[3]
        dt = header[4]
        
        data = np.genfromtxt(str(n) + ".csv", delimiter=",", skip_header=4)
        X = data[:,0:3]
        u = data[:,3:6]
        x = X+u
        triads = data[:,6:15]
        assert N==x.shape[0]

        
        ax.cla()
        ax.plot3D(x[:,0],x[:,1],x[:,2],".-",color="black")
        
        for i in range(N):
            assert triads[i,:].shape == (9,)
            U =  triads[i,:].reshape(3,3)
            plot_triad(ax,x[i,:],U)
            
        
        plt.axis("equal")
        # title = "n={: <6}, t={: <6}, dt={: <6}".format(n, t, dt)
        # plt.title(title)
        plt.title("n="+str(n)+", t="+str(t)+", dt="+str(dt))
        ax.set_xlabel("x")
        ax.set_ylabel("y")
        ax.set_zlabel("z")
        plt.pause(0.1)
    plt.show()
