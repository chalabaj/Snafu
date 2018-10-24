import numpy as np
import os
import shutil
try:
    from errors import error_exit
    from constants import *
except ImportError as ime:
    print("Module {} not found.".format(ime.name))
    exit(1)

def check_restart(restart, cwd):
    rst_file = "restart.in"
    rst_file_path = os.path.join(cwd, rst_file)
    if restart < 0:
            error_exit(11, "(restart < 0)")
    elif restart == 0:
        if (os.path.isfile(rst_file_path)):
            print("File restart.in exists, but restart option is turned off.",
                  "\nRenaming restart.in to restart.in_old.")
            try:
                os.rename("restart.in", "restart.in_old")
            except OsError:
                print("Error, when renaming the restart file. Exiting...")
                exit(1)
        else:
            print("Restart turned off.")
    elif restart == 1:
        rst_file = "restart.in"  
        rst_file_path = os.path.join(cwd, rst_file)
        if (os.path.isfile(rst_file_path)):
            print("Restart={}, start from the last completed".format(restart),
                  "step.\nRestart file {} found.".format(rst_file))
        else:
            error_exit(10, "({})".format(rst_file))
    elif restart > 1:
        rst_file = "restart_{}.in".format(restart)
        rst_file_path = os.path.join(cwd, rst_file)
        if (os.path.isfile(rst_file_path)):
            print("Restart={}, start from {} step".format(restart, restart),
                  "Restart file {} found".format(rst_file))
        else:
            error_exit(10, "({})".format(rst_file))
    return(rst_file_path)
    
def read_restart(rst_file_path):
    #with open restart.in as rsf:
    rst_file = "restart.in"
    with open(rst_file_path, 'r') as rstf:
        step = rstf.readline().split()[1]
        print(step)
    rstf.closed
    pot_eners_array = np.loadtxt(rst_file, dtype=np.float64, delimiter=None, skiprows=8)

    x = np.loadtxt(rst_file, dtype=np.float64, delimiter=None, 
                   skiprows=11, usecols=1)  
    y = np.loadtxt(rst_file, dtype=np.float64, delimiter=None, 
                   skiprows=11, usecols=2)  
    z = np.loadtxt(rst_file, dtype=np.float64, delimiter=None, 
                   skiprows=11, usecols=3)  

    vx = np.loadtxt(rst_file, dtype=np.float64, delimiter=None, 
                   skiprows=15, usecols=1)  
    vy = np.loadtxt(rst_file, dtype=np.float64, delimiter=None, 
                   skiprows=15, usecols=2)  
    vz = np.loadtxt(rst_file, dtype=np.float64, delimiter=None, 
                   skiprows=15, usecols=3)

    fx = np.loadtxt(rst_file, dtype=np.float64, delimiter=None, 
                   skiprows=1�, usecols=1)  
    fy = np.loadtxt(rst_file, dtype=np.float64, delimiter=None, 
                   skiprows=1�, usecols=2)  
    fz = np.loadtxt(rst_file, dtype=np.float64, delimiter=None, 
                   skiprows=1�, usecols=3)             
    print(x)
    print(y)
    print(z)
    print(vx)
    print(vy)
    print(vz)
    print(fx)
    print(fy)
    print(fz)
    
    return()

def print_restart(
        step, time, natoms, at_names, state, timestep,
        x, y, z, vx, vy, vz, fx, fy, fz,
        Ekin, Epot, Etot, Etot_init, pot_eners_array, restart_write):

    inits_line = ("Step: {:d}".format(step),
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
        rsf.write('\n'.join(inits_line))

        np.savetxt(rsf, pot_eners_array, fmt="%20.10f", delimiter=' ', newline='\n')

        rsf.write("Positions AT X Y Z:\n")
        xx = [xi*BOHR_ANG for xi in x.tolist()]
        yy = [yi*BOHR_ANG for yi in y.tolist()]
        zz = [zi*BOHR_ANG for zi in z.tolist()]
        for iat in range(0, natoms):     
            p_line = "{} {:20.10f} {:20.10f} {:20.10f}\n".format(at_names[iat],
                                                                 xx[iat],
                                                                 yy[iat],
                                                                 zz[iat])
            rsf.write(p_line)
            
        rsf.write("Velocities: AT VX VY VZ:\n")
        vvx = vx.tolist()
        vvy = vy.tolist()
        vvz = vz.tolist()
        for iat in range(0, natoms):     
            v_line = "{} {:20.10f} {:20.10f} {:20.10f}\n".format(at_names[iat],
                                                                 vvx[iat],
                                                                 vvy[iat],
                                                                 vvz[iat])
            rsf.write(v_line)

        rsf.write("Forces: AT FX FY FZ:\n")
        ffx = fx.tolist()
        ffy = fy.tolist()
        ffz = fz.tolist()
        for iat in range(0, natoms):
            f_line = "{} {:20.10f} {:20.10f} {:20.10f}\n".format(at_names[iat],
                                                               ffx[iat],
                                                               ffy[iat],
                                                               ffz[iat])
            rsf.write(f_line)
    rsf.closed

    if not (step%restart_write):
        rst_file = "restart_{}.in".format(step)
        print("Writing restart information to {} file.".format(rst_file))
        shutil.copy("restart.in", rst_file, follow_symlinks=True) 

    #print(pot_eners_array)
    #for irow in range(0, np.size(pot_eners_array, 0)):
        #line2 = ("Epot_array: \n",
        #         "{5.8f}".format(' '.join(str(pot_eners_array[irow,:]))))
    
    #print(np.array2string(pot_eners_array,precision=10))
          #formatter={'float_kind':lambda x: "{}".format(x)}))
      # line2 = ("".join("{2s} {5.8f} {5.8f} {5.8f}\n".format(at_names[iat],vx[iat] ,vy[iat],vz[iat])))        
    return()
    
