"""
SNAFU code
SNAFU can do ab initio nonadiabatic MD withot calculating non-adiabatic couplings
Two scheme with/without diabatization

1) adiabatic-potential-based formula derived within Landau-Zener (APLZ):
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
from snafu.init   import read_geoms, read_velocs
from snafu.init   import create_output_file, init_forces
from snafu.masses import assign_masses
from snafu.errors import error_exit
from snafu.propagate import update_velocities, update_positions
from snafu.print import print_positions
 
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
#----------------INIT-------------------------------------------------------------------------------

if __name__ == "__main__":
    
    print("Starting SNAFU.\n")
    print("Python exec: ",sys.base_exec_prefix,"ver:",sys.version[:5])
    print("System platform: ",sys.platform)  
    print("Working directory: ",cwd)
    print(liner)
    
# FILE CHECK - OBTAIN PATH TO FILES
    input_file_path, geom_file_path, veloc_file_path, veloc_init = file_check(cwd)
    
#READ INPUT VARIABLE (read as strings) - SET THEM AS GLOBAL VARIABLES: 
    input_vars = read_input(input_file_path)
    globals().update(input_vars)  
    natoms   = int(natoms)                                     # loaded variable are all strings
    maxsteps = int(maxsteps)
    state = int(init_state)                                    # Current state during propagation (for t =0 => init state)
    dt = float(timestep)
    
#READ INITIAL GEOMETRY AND VELOCITIES AND CREATE ARRAYS FOR FORCES:
    # each array (x,y,z) hold NATOMS positions, 
    # e.g. x[0] is x position of 1st atom, fx[1] force in x direction on 2 atom, etc.
    at_names,x,y,z = read_geoms(natoms,geom_file_path)         # names of atomos, x,y,z positions (1D)
    vx,vy,vz = read_velocs(veloc_init,natoms,veloc_file_path)  # velocity array (1D) 
    fx,fy,fz,fx_new,fy_new,fz_new = init_forces(natoms)         # Create empty forces array (1D)
    
#OBTAIN MASSES:
    masses = assign_masses(at_names) 
    am = [ mm * amu for mm in masses]                          # atomic mass units conversion
    print(liner)
    if debug == 1: print("atomic masses:\n",am)
    print("Molecular systems:\nAt  Mass     X     Y     Z:")
    for iat in range(0,natoms):
        print("".join("%2s" " " "%2.2f"  %(at_names[iat], masses[iat]))," %2.4f %2.4f %2.4f"  %(x[iat],y[iat],z[iat]))  # just nice output print
    print(liner)     
    
#CREATE OUTPUT FILES:
    # where to store propagated position, velocities, observables
    files = [ "energies.dat", "velocities.dat", "gradients.dat", "movie.xyz", "restart.dat" ]  
   
    create_output_file(files)

#---------------INIT DONE-------------------------------------------------------------------    
# CALC INITIAL ENERGIES AND GRADIENTS
    # fx,fy,fz = calc_force(state,) # position at current step
 
     
# MAIN LOOP 
    #center of mass reduction TODO
    for step in range(1,maxsteps):
        
        if debug == 1: print("fx:\n ",fx)
        
        x, y, z = update_positions(natoms,dt,am,x,y,z,vx,vy,vz,fx,fy,fz)   # new positions (t+dt)
        
        # fx_new, fy_new, fz_new = calc_force(x,y,z,state)                 # calc forces for new positions
        
        vx, vy, vz = update_velocities(natoms,dt,am,vx,vy,vz,fx,fy,fz,fx_new,fy_new,fz_new) # propagate velocities using new forces
        
        fx = fx_new
        fy = fy_new
        fz = fz_new
        
        #alcc_hop
        # vel_adjustment
        #calc_energies
        #print_info(step, pos, ener)
        time = step * au_fs 
        # print_positions(step,time,natoms, at_names, x, y, z):   
    

    #prepare files - energies, vel, xyz pos for production data, if exists and rstart = 0 then crash.


