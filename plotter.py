import numpy as np
import matplotlib.pyplot as plt
import os
from mpl_toolkits.mplot3d import Axes3D

# fmt: off

#need this to make ctrl+c close plot window
import signal
signal.signal(signal.SIGINT, signal.SIG_DFL)

SMALL_VAL = 1e-6
MAX_TRIADS = 10
MAX_NODES=1000

#These are 1-indexed:
LINE_HEADER_DATA = 2
LINE_ENERGY_DATA = 4
FIRST_LINE_OUPUT_DATA=6



def plot_triad(ax, x,U,scale=1):
    assert x.shape==(3,)
    assert U.shape==(3,3)
    U*=scale
    ax.quiver(x[0],x[1],x[2],U[0,0],U[1,0],U[2,0],color="r")
    ax.quiver(x[0],x[1],x[2],U[0,1],U[1,1],U[2,1],color="g")
    ax.quiver(x[0],x[1],x[2],U[0,2],U[1,2],U[2,2],color="b")

def plot_hole(x,r,ax):
    color="grey"
    N_hole = x.shape[0]
    theta = np.linspace(0,2*np.pi,20)
    e2 = np.array([0.,1.,0.])
    for i in range(N_hole-1):
        Rx = r[i]*np.cos(theta)
        Ry = r[i]*np.sin(theta)
        Z = np.ones_like(theta)*x[i,0]
        ax.plot3D(Z,Rx,Ry,color=color)
        x1 = np.array([x[i,0],x[i+1,0]])
        x2 = np.array([x[i,1],x[i+1,1]])+r[i]
        x3 = np.array([x[i,2],x[i+1,2]])
        ax.plot3D(x1,x2,x3,color=color)

component_to_index = {"x":0,"y":1,"z":2,"all":np.arange(0,3)}
vec3_data_first_index = {"u":3, "v": 15, "omega_u": 18} #add more when needed

def get_vector_components_from_data(data,datatype,component="all"):
    assert datatype in vec3_data_first_index
    assert component in component_to_index
    first_index = vec3_data_first_index[datatype]
    return data[:,first_index+component_to_index[component]]
   

class Plotter:    
    def __init__(self,write_gif=True) -> None:
        os.chdir(os.getcwd())
        print("Python script running from the following directory:\n",os.getcwd())
        self.output_dir = os.path.join(os.getcwd(),"output")
        try:
            header = self.read_header(0)
        except:
            print("No output in output directory \"" + self.output_dir + "\", remember to set save_csv=true")
            exit(1)
        self.N = int(header[0])
        data = self.read_data(0)
        X = data[:,0:3]
        self.L0 =  X[-1][0] - X[0][0]
        self.n_steps = int(header[1])
        self.n_write = int(header[2])
        self.n_plot_triad = max(int(self.N/MAX_TRIADS),1)
        self.check_energy_balance = bool(int(header[5]))
        self.contact_enabled = bool(int(header[6]))
        self.write_gif=write_gif
        if write_gif:
            self.output_tmp = "output-tmp"
            if os.path.exists(self.output_tmp) is False:
                os.makedirs(self.output_tmp)
        if self.N>MAX_NODES:
            node_stride = int(self.N/MAX_NODES)
            self.i_nodes = np.arange(0,self.N,node_stride)
        else:
            self.i_nodes = np.arange(0,self.N,1)  
        
        self.read_borehole()
        
    def read_header(self,n):
        return np.genfromtxt(os.path.join(self.output_dir, str(n)+".csv"), delimiter=",", skip_header=LINE_HEADER_DATA-1,max_rows=1)
    def read_data(self,n):
        data =  np.genfromtxt(os.path.join(self.output_dir, str(n) + ".csv"), delimiter=",", skip_header=FIRST_LINE_OUPUT_DATA-1)
        return data
    def read_energy_balance(self,n):
        assert self.check_energy_balance
        data = np.genfromtxt(os.path.join(self.output_dir, str(n)+".csv"), delimiter=",", skip_header=LINE_ENERGY_DATA-1,max_rows=1)
        assert data.shape == (4,)
        KE = data[0]
        W_int = data[1]
        W_ext = data[2]
        E_res = data[3]
        return KE, W_int,W_ext,E_res
            
    def read_header_transient(self,n):
        header = self.read_header(n)
        self.t = header[3]
        self.dt = header[4]
        
    def read_borehole(self):
        if self.contact_enabled:
            file=os.path.join(self.output_dir,"borehole.csv")
            #header
            data = np.genfromtxt(file,delimiter=",",skip_header=1,max_rows=1)
            self.N_hole=data
            #hole geometry
            data = np.genfromtxt(file,delimiter=",",skip_header=3)
            assert data.shape[0]==(self.N_hole)
            self.x_hole = data[:,0:3]
            self.r_hole = data[:-1,3]
        
        
    def animate_3d(self):
        ax = plt.figure(figsize=(8,6)).add_subplot(111,projection="3d")
        for n in range(self.n_steps):
            if n % self.n_write != 0: continue
            self.read_header_transient(n)
            
            data = self.read_data(n)
            X = data[:,0:3]
            
            quiver_scale = self.L0/10
            u = data[:,3:6]
            x = X+u
            assert self.N==x.shape[0]
            triads = data[:,6:15]
            100000000
            ax.cla()
            ax.plot3D(x[self.i_nodes,0],x[self.i_nodes,1],x[self.i_nodes,2],".-",color="black")
            
            for i in range(self.N):
                if i%self.n_plot_triad!=0: continue
                assert triads[i,:].shape == (9,)
                U =  triads[i,:].reshape(3,3)
                plot_triad(ax,x[i,:],U,quiver_scale)
                
            plt.axis("equal")
            # title = "n={: <6}, t={: <6}, dt={: <6}".format(n, t, dt)
            # plt.title(title)
            if self.contact_enabled: plot_hole(self.x_hole,self.r_hole,ax)
           
            plt.title("n = %d, t = %.3f [s], dt = %.6f [s]" % (n,self.t,self.dt))
            
            
            axlen = 1.2*self.L0
            #axlen = 0.7*self.L0
            ax.set_xlim(0,axlen)
            ax.set_ylim(-axlen/2,axlen/2)
            ax.set_zlim(-axlen/2,axlen/2)
            ax.set_xlabel("x")
            ax.set_ylabel("y")
            ax.set_zlabel("z")
            ax.set_xlim()   
            if self.write_gif:
                filepath = os.path.join(self.output_tmp, str(n) + ".png")
                plt.savefig(filepath,dpi=200)
            plt.pause(SMALL_VAL)
            
        if self.write_gif: self.create_gif()
    
    def get_end_node_disp_transient(self):
        
        n_vals = (self.n_steps - 1) // self.n_write + 1
      
        t = np.zeros((n_vals,))
        u_tip = np.zeros((n_vals,3))
        i = 0
        for n in range(0,self.n_steps,self.n_write):
            if n % self.n_write != 0: continue
            self.read_header_transient(n)
            t[i] = self.t
            
            data = self.read_data(n)
            disp_all_nodes = get_vector_components_from_data(data,"u")
            u_tip[i,:] = disp_all_nodes[-1,:]
            i += 1
        assert n_vals==i
        return t,u_tip
            
    
    def plot_end_node_transient(self,components):
        
        assert len(components)==len(set(components))
        for c in components:
            assert c=="x" or c=="y" or c=="z" 
            
        t,u_tip = self.get_end_node_disp_transient()

        for i,c in enumerate(components):
            plt.figure()
            print("Max |"+c+"-displacement| = "+str(max(abs(u_tip[:,i]))))
            plt.plot(t,u_tip[:,i])
            plt.grid()
            plt.tight_layout()
            plt.xlabel("t [s]")
            plt.ylabel(c+"-displacement [m]")
            
    def animate_vertical_disp(self):
        plt.figure()
        for n in range(self.n_steps):
            if n % self.n_write != 0: continue
            self.read_header_transient(n)
            data = self.read_data(n)
            u3 = data[:,5]
            plt.cla()
            plt.plot(u3,'-')
            plt.xlabel("i")
            plt.ylabel("u3")
            plt.ylim([-self.L0/5,self.L0/5])
            plt.pause(SMALL_VAL)

    def plot_specific_kinetic_energy_component_wise(self):
        steps = np.arange(0,self.n_steps,self.n_write)
        t = np.zeros_like(steps, dtype=float)
        KE_t = np.zeros_like(t)
        KE_w_u = np.zeros_like(t)
        KE_w_u_1 = np.zeros_like(t)
        KE_w_u_2 = np.zeros_like(t)
        KE_w_u_3 = np.zeros_like(t)
        for i,n in enumerate(steps):
            self.read_header_transient(n)
            t[i] = self.t
            data = self.read_data(n)
            assert data.shape[1] == 21
            v1 = data[:,15]
            v2 = data[:,16]
            v3 = data[:,17]
            w_u_1 = data[:,18]
            w_u_2 = data[:,19]
            w_u_3 = data[:,20]
            KE_t[i] = 0.5*np.sum(v1**2 + v2**2 + v3**2)
            KE_w_u_1[i] = 0.5*np.sum(w_u_1**2)
            KE_w_u_2[i] = 0.5*np.sum(w_u_2**2)
            KE_w_u_3[i] = 0.5*np.sum(w_u_3**2)
            KE_w_u[i] = KE_w_u_1[i] + KE_w_u_2[i] + KE_w_u_3[i]
        plt.figure()
        # plt.plot(t, KE_t, label="KE_t")
        # plt.plot(t, KE_w_u, label="KE_w_u")
        plt.plot(t, KE_w_u_1, label="KE_w_u1")
        plt.plot(t, KE_w_u_2, label="KE_w_u2")
        plt.plot(t, KE_w_u_3, label="KE_w_u3")
        plt.ylim([0, 10000])
        plt.legend()
        plt.xlabel("t[s]")
        plt.ylabel("Specific kinetic energy")
        
    def animate_omega_u(self):
        plt.figure()
        steps = np.arange(0,self.n_steps,self.n_write)
        i = np.arange(0,self.N,1)
        for n in steps:
            self.read_header_transient(n)
            data = self.read_data(n)
            w_u_1 = data[:,18]
            w_u_2 = data[:,19]
            w_u_3 = data[:,20]
            #v_1 = data[:,15]
            plt.cla()
            plt.plot(w_u_1,'.-',label="omega_u_1")
            plt.plot(w_u_2,'.-',label="omega_u_2")
            plt.plot(w_u_3,'.-',label="omega_u_3")      
            #plt.plot(v_1,'.-')
            plt.xlabel("i")
            plt.ylabel("omega_u [rad/s]")
            plt.title("t="+str(self.t))
            limval = 1
            limval=1000
            plt.ylim([-limval,limval])
            plt.tight_layout()
            #plt.ylim([-self.L0/5,self.L0/5])
            plt.pause(SMALL_VAL)

    def plot_energy_balance(self):
        if not self.check_energy_balance:
            print("Energy balance was not monitored")
            return
        
        t = []
        KE = []
        W_int = []
        W_ext = []
        E_res = []
        for n in range(self.n_steps):
            if n % self.n_write != 0: continue
            self.read_header_transient(n)
            t.append(self.t)
            KE_, W_int_,W_ext_, E_res_ = self.read_energy_balance(n)
            KE.append(KE_)
            W_int.append(W_int_)
            W_ext.append(W_ext_)
            E_res.append(E_res_)
        plt.figure()
        t = np.array(t)
        KE = np.array(KE)
        W_int = np.array(W_int)
        W_ext = np.array(W_ext)
        E_res = np.array(E_res)
        plt.plot(t,KE,label=r"$K$")
        plt.plot(t,W_int,label=r"$W^{int}$")
        plt.plot(t,W_ext,label=r"$W^{ext}$")
        plt.plot(t,E_res,label=r"$E^{res}$")
        #lim=20000
        #plt.ylim(0,lim)
        #plt.ylim(-lim/2,lim/2)
        plt.legend()
        plt.xlabel("$t$ [s]")
        plt.ylabel("Energy $[J]$")
        
        
    def create_gif(self):
        assert self.write_gif
        print("Writing Gif...")
        from PIL import Image

        image_files = [str(n) + ".png" for n in range(0, self.n_steps, self.n_write)]
        images = []
        for img_file in image_files:
            img = Image.open(os.path.join(self.output_tmp, img_file))
            images.append(img)

        output_name = "beam.gif"
        dt_frame = 1000.0 / 60  # 60 frames per second
        images[0].save(
            output_name, save_all=True, append_images=images[1:], duration=dt_frame, loop=0
        )
        for img in images:
            img.close()
        import shutil
        shutil.rmtree(self.output_tmp)
        print("Gif done")
            

if __name__ == "__main__":
    
    p = Plotter(write_gif=False)
    #p.plot_specific_kinetic_energy_component_wise()
    #p.plot_end_node_transient()
    p.plot_energy_balance()
    #p.animate_vertical_disp()  
    p.animate_3d()
    #p.animate_omega_u()
    plt.show()
