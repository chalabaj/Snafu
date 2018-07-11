"""
SNAFU code
SNAFU can do ab initio nonadiabatic MD withot calculating non-adiabatic couplings
Two scheme with/without diabatization

1) adiabatic-potential-based formula derived within Landau-Zener:
AK Belyaev, PHYSICAL REVIEW A 84, 014701 (2011),DOI: 10.1103/PhysRevA.84.014701
AK Belyaev, The Journal of Chemical Physics 147, 234301 (2017); doi: 10.1063/1.5000718

2) diabatization scheme: Le Yu, Phys.Chem.Chem.Phys., 2014, 16, 25883
    
TO DO: step back if hop, restart file will contain geometry, velocities, curent state, timestep atd...
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
from snafu.init   import read_geoms, read_velocs,create_output_file 
from snafu.masses import assign_masses
from snafu.errors import error_exit
from snafu.propagate import verlet_step
 
# Constants
au_fs = 0.02418884326505      #atomic units to femtosecs
au_eV = 27.21139
amu   = 1822.8885             # atomic mass unit  me = 1 AMU*atomic weight
ang_bohr = 1.889726132873     # agstroms to bohrs

# Observed variable
aEk = 0.0   # average kinetic energy
cPE = 0.0   # current potential energy to save
cKE = 0.0   # current kinetic energy to save

liner = ("_______________________________________________")
debug = 1   
if __name__ == "__main__":
    
    print("Starting SNAFU.\n")
    print("Python exec: ",sys.base_exec_prefix,"ver:",sys.version[:5])
    print("System platform: ",sys.platform)  
    print("Working directory: ",cwd)
    print(liner)
    
    # FILE CHECK - OBTAIN PATH TO FILES
    input_file_path, geom_file_path, veloc_file_path, veloc_init = file_check(cwd)
    
    #READ INPUT VARIABLE (read as strings) - SET THEM AS GLOBAL VARIABLES 
    input_vars = read_input(input_file_path)
    globals().update(input_vars)  
    natoms   = int(natoms)                                     # load variable are all strings
    maxsteps = int(maxsteps)
    state = int(init_state)                                    # current state during propagation (t =0 => init state)
    dt = float(timestep)
    
    #READ INITIAL GEOMETRY AND VELOCITIES
    at_names,x,y,z = read_geoms(natoms,geom_file_path)         # names of atomos, x,y,z position - vectors x[0] is x position of the 1st atom
    vx,vy,vz = read_velocs(veloc_init,natoms,veloc_file_path)  # vectors vx[0] is x velocity of the 1st atom
    
    #OBTAIN MASSES
    masses = assign_masses(at_names) 
    am = [ mm * amu for mm in masses]                          # atomic mass units conversion
    print(liner)
    if debug == 1: print("atomic masses:\n",am)
    print("Molecular systems:\nAt  Mass     X     Y     Z:")
    for iat in range(0,natoms):
        print("".join("%2s" " " "%2.3f"  %(at_names[iat], masses[iat]))," %2.3f %2.3f %2.3f"  %(x[iat],y[iat],z[iat]))  # just nice output print
    print(liner)     
    
    #CREATE OUTPUT FILES
    files = [ "energies.dat", "velocities.dat", "gradients.dat", "movie.xyz", "restart.dat" ]  # where to store propagate position, velocities
    create_output_file(files)
 
    # MAIN LOOP
    #center of mass reduction TODO:
    for step in range(1,maxsteps):
        fx = np.zeros(shape=(natoms,1))+1  # TO DO CALL FORCES
        if debug == 1: print("fx:\n",fx)
        #velocity verlet step
        x_new = verlet_step(natoms,dt,am,x,vx,fx)
        #alcc_hop
        # vel_adjustment
        #calc_energies
        #print_info(step, pos, ener)
        time = step * au_fs 
        # store data   
    
def print_pos():
    line = "  ".join("%1d" "%1d" "%1d" %(x[1], y[1], z[1]))
    print(line)
    return
    #prepare files - energies, vel, xyz pos for production data, if exists and rstart = 0 then crash.


