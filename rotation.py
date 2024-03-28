import numpy as np
import matplotlib.pyplot as plt

import signal

signal.signal(signal.SIGINT, signal.SIG_DFL)

# fmt: off

SMALL_SCALAR = 0.3



def is_orthogonal(R):
    diff = np.linalg.inv(R) - R.T 
    for i in range(3):
        for j in range(3):
            if (abs(diff[i,j])>SMALL_SCALAR):
                print("R^-1 - R_T =",diff)
                return False
    if (abs(np.linalg.det(R)-1)>SMALL_SCALAR):
        print("det R =",np.linalg.det(R))
        return False
    else: return True
    

def skew(a):
    S = np.array([[0, -a[2], a[1]],
                [a[2], 0, -a[0]],
                [-a[1], a[0], 0]])
    return S



n_max=1000



I = np.identity(3)
U = np.identity(3)
Ju = np.array([0.322136, 1718.06, 1718.06])
#Ju = np.array([1.,1000.,1000.])
omega = np.array([100., 1., 0])
omega_u = U.T@omega
omega_u_dot = np.array([0.,0.,0.])

dt = 0.00168375

KE = np.zeros((n_max,))
times = np.arange(0.,dt*n_max,dt)
assert np.size(times)==n_max
dt /=10

if __name__=="__main__":
    
    
    print("omega at t0 = ",omega, ", norm = ",np.linalg.norm(omega))
    
    
    for n in range(n_max):   
        t = n*dt
        kinetic_energy = 0.5*(Ju[0]*omega_u[0]**2 + Ju[1]*omega_u[1]**2 + Ju[2]*omega_u[2]**2)
        print("KE at n=",n," =", kinetic_energy) 
        KE[n] = kinetic_energy
        M = np.array([0., 0., 0.])
        omega_u_dot[0] = (M[0] - (Ju[2]-Ju[1])*omega_u[1]*omega_u[2])/Ju[0]
        omega_u_dot[1] = (M[1] - (Ju[0]-Ju[2])*omega_u[0]*omega_u[2])/Ju[1]
        omega_u_dot[2] = (M[2] - (Ju[1]-Ju[0])*omega_u[0]*omega_u[1])/Ju[2]
        
        omega_u += dt*omega_u_dot
        
        theta = np.linalg.norm(omega_u*dt)
        print("theta [deg] =",theta*180/np.pi)
        
        
        U = U*np.linalg.inv(I - 0.5*dt*skew(omega_u))@(I + 0.5*dt*skew(omega_u))
        assert is_orthogonal(U)
        
        print("omega at t = ",t," = ",omega, ", norm = ",np.linalg.norm(omega))

    plt.plot(times,KE)
    plt.xlabel("t [s]")
    plt.ylabel("Kinetic energy [J]")
    plt.ylim([0., np.max(KE)*1.2])
    plt.show()
