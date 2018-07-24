"""
Main code. Calls: initialization routines (geometry, velocities)
calls for files check
calls for array init
main iteration loop over stepsr
"""

# SYSTEM IMPORTS:
import math
import sys
import os
import random
import time
import re
import numpy as np
np.set_printoptions(precision=8) # needed for accuracy and stability testing only

from datetime import datetime

# ENVIRONMENT LAYER to find local modules
# .bashrc should contain PYTHONPATH=/path/to/SNAFU/snafu
# or sys.path.append("/path/to/SNAFU")
# or export PYTHONPATH=/path/to/SNAFU/snafu before code calls

cwd = os.getcwd()
#sys.path.append(cwd)
#sys.path.append('/home/chalabaj/SNAFU/snafu')
#sys.path.append('/home/chalabaj/SNAFU')

# LOCAL IMPORT
modules_files = [
    'init.py',
    'masses.py',
    'errors.py',
    'prints.py',
    'propagate.py',
]
def load_local_modules():
    from init import file_check, read_input
    from init import read_geoms, read_velocs
    from init import create_output_file, init_forces_potenergs
    from init import com_removal
    from masses import assign_masses
    from errors import error_exit
    from prints import print_positions, print_velocities
    from prints import print_snafu
    from propagate import update_velocities, update_positions
    from propagate import calc_forces, calc_energies
    return()
    
try:
    env_set = True
    env_set = False for pathline in sys.path for mod_file in modules_files if os.path.isfile(os.path.join(pathline,
    
    load_local_modules() 
except:    
else:
    print("All modules loaded succesfully.")


# Constants
au_fs = 0.02418884326505e0      #atomic units to femtosecs
au_eV = 27.21139
amu   = 1822.888484264545e0             # atomic mass unit  me = 1 AMU*atomic weight
ang_bohr = 1.889726132873e0     # agstroms to bohrs


# Observed variable
step = 0
dE = 0.0   # energy change since initial energies
Etot_init = 0.0  # setting variable , total energy at the beginning 
time = 0.0
liner = ("_")*70
debug = 1   
#----------------INIT-------------------------------------------------------------------------------

if __name__ == "__main__":

    print("Starting:\n")
    print_snafu()
    startTime = datetime.now()
    print("Simulation started at: {}".format(startTime))
    print("Python exec: ",sys.base_exec_prefix,"ver:",sys.version[:5])
    print("System platform: ",sys.platform)  
    print("Working directory: ",cwd)
    print(liner)
    
# FILE CHECK - OBTAIN PATH TO FILES
    input_file_path, geom_file_path, veloc_file_path, veloc_init = file_check(cwd)
    
#READ INPUT VARIABLE (read as strings) - SET THEM AS GLOBAL VARIABLES: 
    input_vars, ab_initio_file_path = read_input(cwd,input_file_path)
    globals().update(input_vars)  
    natoms   = int(natoms)                                     # loaded variable are all strings
    maxsteps = int(maxsteps)
    state = int(init_state)                                    # Current state during propagation (for t =0 => init state)
    dt = float(timestep)
    nstates = int(nstates)
    
#READ INITIAL GEOMETRY AND VELOCITIES AND CREATE ARRAYS FOR FORCES:
    # each array (x,y,z) hold NATOMS positions, 
    # e.g. x[0] is x position of 1st atom, fx[1] force in x direction on 2 atom, etc.
    at_names, x, y, z = read_geoms(natoms,geom_file_path)                           # names of atomos, x,y,z positions (1D)
    vx, vy, vz = read_velocs(veloc_init,natoms,veloc_file_path)                    # velocity array (1D) 
    fx, fy, fz, fx_new, fy_new,fz_new, pot_eners = init_forces_potenergs(natoms,nstates)         # Create empty forces and potential energy arrays (1D)
    
#OBTAIN MASSES:
    masses = assign_masses(at_names) 
    am = [ mm *  amu for mm in masses]                          # atomic mass units conversion
    print(liner)
    #if debug == 1: print("atomic masses:\n",am)
    #print("Molecular systems:\nAt  Mass     X     Y     Z:")
    for iat in range(0,natoms):
        print("".join("%3s" " " "%3.3f"  %(at_names[iat], masses[iat]))," %3.6f %2.6f %2.6f"  %(x[iat],y[iat],z[iat]))  # just nice output print
    print(liner)     
    
#CREATE OUTPUT FILES:
    # where to store propagated position, velocities, observables
    # geom.dat hold current geometry for which to compute E, grads  
    #files = [ "energies.dat", "velocities.xyz", "gradients.dat", "movie.xyz" ]  
    #create_output_file(files)

#---------------INIT DONE-------------------------------------------------------------------    

# CENTER OF MASS REMOVAL
    #x, y, z = com_removal(x,y,z,am)
    
# CALC INITIAL ENERGIES AND GRADIENTS
    print("Step      Time/fs     Energy change from start/eV  Hoppping") 
    
    fx, fy, fz, pot_eners = calc_forces(step, natoms, at_names, state, nstates, ab_initio_file_path, x, y, z, fx, fy, fz, pot_eners) # position at current step   
    
    Ekin, Epot, Etot, dE = calc_energies(step, time, natoms, am, state, pot_eners, vx, vy, vz, Etot_init) #Etot_init = 0
    Etot_init = Etot  # Total energy at the beginning to calc. energy changes during propagation
    
    
# MAIN LOOP 
    #center of mass reduction TODO
    for step in range(1,maxsteps+1):
        
        x, y, z = update_positions(natoms,dt,am,x,y,z,vx,vy,vz,fx,fy,fz)                                                                  # new positions (t+dt)
       
        fx_new, fy_new, fz_new, pot_eners = calc_forces(step,natoms, at_names, state, nstates, ab_initio_file_path, x, y, z, fx_new, fy_new, fz_new, pot_eners)  # calc forces for new positions
       
        vx, vy, vz = update_velocities(natoms,dt,am,vx,vy,vz,fx,fy,fz,fx_new,fy_new,fz_new) # propagate velocities using new forces
        
        fx = np.copy(fx_new) # copied list instead of just referencing
        fy = np.copy(fy_new)
        fz = np.copy(fz_new)
        
        #if hopping == "1":
          #alcc_hop
        # vel_adjustment
        
        hop = "No"
        time = step * dt*  au_fs 
        Ekin, Epot, Etot, dE = calc_energies(step, time, natoms, am, state, pot_eners, vx, vy, vz, Etot_init)
        
        print(" {:<3d} {:>10.2f} {:>20.4e} {:>20s}\n".format(step,time,dE* au_eV, hop))
        
     # save positions and velocities to movie and velocity file
        print_positions(step,time,natoms, at_names, x, y, z)   
        print_velocities(step,time,natoms, at_names, vx, vy, vz)
      
    
    print("JOB completed.") 
    print(liner)
    
    stopTime = datetime.now()
    simtime = (datetime.now() - startTime)
    print("Simulation ended at: {}".format(stopTime))
    print("Overall simulation time (hh:mm:ss): {}".format(simtime))     

