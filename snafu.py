"""
SNAFU code
SNAFU can do ab initio nonadiabatic MD withot calculating non-adiabatic couplings
Two scheme with/without diabatization

1) adiabatic-potential-based formula derived within Landau-Zener:
AK Belyaev, PHYSICAL REVIEW A 84, 014701 (2011),DOI: 10.1103/PhysRevA.84.014701
AK Belyaev, The Journal of Chemical Physics 147, 234301 (2017); doi: 10.1063/1.5000718

2) diabatization scheme: Le Yu, Phys.Chem.Chem.Phys., 2014, 16, 25883
"""
import math
import sys, os
import numpy as np
import random
import time

cwd = os.getcwd()
sys.path.append(cwd)
from check import file_check, read_input, error_exit

# Constants
au_fs = 0.02
au_eV = 27.21139

# Observed variable
aEk = 0.0   # average kinetic energy
cPE = 0.0   # current potential energy to save
cKE = 0.0   # current kinetic energy to save



def main_loop():
    """
    Main rutine to propagate atoms
    Calling external code to calc. force, energies
    Calcul transition probabilities 
    Readjust velocities
    Check. energy conservation
    """
#  x_old  x(t-dt)
#  x_     x(t)
#  x_new  x(t+dt)
    #cacl_force
    for step in range (o,maxsteps):
        print("verlet")
        #verlet_step
        #calc_hop
        # vel_adjustment
        #calc_energies
        #print_info(step, pos, ener)
        
             
        
    return

def verlet_step(x,y,z,m,t):
    for iat in range(1,natoms+1):   # upper index excluded
     print("for atom")   
    #eq 1 px(iat,t + 1/2 dt) = px(iat,t) + 1/2*dt*fi(iatt))
    #eq 2
    #call forces
    # eq 3
    return()

def read_geoms_veloc(geom_file_path,veloc_file_path):
    #if restart == 1 : read  last two geoms
    # could be xyz matrix [natoms,3] as well, but this should be more clear reading
    x = np.zeros(shape=(natoms,1)) 
    y = np.zeros(shape=(natoms,1)) 
    z = np.zeros(shape=(natoms,1)) 
 # READ INITIAL POSITIONS:   
    with open(geom_file_path,'r') as igf:  # igf input geom file 
     atoms = igf.readline()  # first line in geom file is number of atoms
     if not (int(atoms) == natoms):  error_exit(2)                         
     garbage = igf.readline()  # comment 
     for iat in range(0,natoms)
        line = igf.readline().split()
        x[iat] = float(line[1])
        y[iat] = float(line[2])
        z[iat] = float(line[3])
     igf.close()
 # READ INITIAL VELOCITIES:
   if
    vx = np.zeros(shape=(natoms,1)) 
    vy = np.zeros(shape=(natoms,1)) 
    vz = np.zeros(shape=(natoms,1)) 
    with open(veloc_file_path,'r') as ivf:  # ivf input veloc file 
     atoms = ivf.readline()  # first line in geom file is number of atoms
     if not (int(atoms) == natoms):  error_exit(2)                         
     garbage = ivf.readline()  # comment 
     for iat in range(0,natoms)
        line = ivf.readline().split()
        vx[iat] = float(line[1])
        vy[iat] = float(line[2])
        vz[iat] = float(line[3])
     ivf.close()
    return
#return x_new, y_new, z_new, v_new, v_new, v_new,
    
if __name__ == "__main__":

    print("Starting SNAFU.\n")
    # File check
    input_file_path, geom_file_path, veloc_file_path, veloc_init = file_check(cwd)
    # Read input parametrs and set them as variable  - all strings => must convert later
    input_vars = read_input(input_file_path)
    globals().update(input_vars)  #Make vars availible globally
    natoms = int(natoms)
    #read initial/restart geometry and velocities
    x,y,z,vx,vy,vz = read_geoms_veloc(geom_file_path,veloc_file_path)    
   
    
    #prepare files - energies, vel, xyz pos for production data, if exists and rstart = 0 then crash.
    # create array  - pos, velocities, energies, gradients
    # load initial pos, velocities
