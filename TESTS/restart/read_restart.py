import numpy as np

def read_restart():
    #with open restart.in as rsf:
    rst_file = "restart.in"
    with open(rst_file, 'r') as rstf:
        step = rstf.readline().split()[1]
        print(step)
    rstf.closed
    a = np.loadtxt(rst_file, dtype=np.float64, delimiter=None, skiprows=8)  
    print(a)
    return()

def print_restart(
        step, time, natoms, at_names, state, timestep,
        x, y, z, vx, vy, vz, fx, fy, fz,
        Ekin, Epot, Etot, Etot_init, pot_eners_array, restart_write):

    infoline = ("Step: {:d}".format(step),
                "State: {:d}".format(state),
                "Natoms: {:d}".format(natoms),
                "Ekin: {:14.10f}".format(Ekin),
                "Epot: {:14.10f}".format(Epot),
                "Etot: {:14.10f}".format(Etot),
                "Etot_init: {:14.10f}".format(Etot_init),
                "Pot_eners_array:\n"
                )

    rst_file = "restart.in"
        
    with open(rst_file, "w") as rsf: 
        rsf.write('\n'.join(infoline))
        np.savetxt(rsf, pot_eners_array, fmt="%20.10f", delimiter=' ', newline='\n')
    rsf.closed
        
    if not (step%restart_write):
        rst_file = "restart_{}.in".format(step)
        print("Writing restart information to {} file.".format(rst_file))
        with open(rst_file, "w") as rsf: 
            rsf.write('\n'.join(infoline))
            np.savetxt(rsf, pot_eners_array, fmt="%20.10f", delimiter=' ', newline='\n')
        rsf.closed 
    #print(pot_eners_array)
    #for irow in range(0, np.size(pot_eners_array, 0)):
        #line2 = ("Epot_array: \n",
        #         "{5.8f}".format(' '.join(str(pot_eners_array[irow,:]))))
    
    #print(np.array2string(pot_eners_array,precision=10))
          #formatter={'float_kind':lambda x: "{}".format(x)}))
      # line2 = ("".join("{2s} {5.8f} {5.8f} {5.8f}\n".format(at_names[iat],vx[iat] ,vy[iat],vz[iat])))        
    return()
    
if __name__ == "__main__":
    read_restart()
    
    exit(0)
    
