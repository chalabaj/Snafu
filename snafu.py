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

#env layer on cluster to find module sys.path.append? /home/XXX/Snafu/snafu
cwd = os.getcwd()
sys.path.append(cwd)

from snafu.init   import file_check, read_input
from snafu.init   import read_geoms, read_velocs 
#from snafu.masses import assign_masses
from snafu.errors import error_exit
 
# Constants
au_fs = 0.02418884326505      #atomic units to femtosecs
au_eV = 27.21139
amu   = 1822.8885             # atomic mass unit  me = 1 AMU*atomic weight
ang_bohr = 1.889726132873     # agstroms to bohrs


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
    at_names,x,y,z = read_geoms(natoms,geom_file_path) 
    vx,vy,vz = read_velocs(veloc_init,natoms,veloc_file_path)
    mass = assign_masses(at_names) 
    mass = 
    (m * amu for m in mass)
    
    print(at_names)
    #masses = assign_mass(names,natoms)
  
def print_pos():
    line = "  ".join("%1d" "%1d" "%1d" %(x[1], y[1], z[1]))
    print(line)
    return
    #prepare files - energies, vel, xyz pos for production data, if exists and rstart = 0 then crash.


