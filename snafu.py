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

# ENVIRONMENT LAYER & LOCAL IMPORT 
# find local modules
# .bashrc should contain SNAFU_DIR=/path/to/SNAFU/snafu
# or launcher exports SNAFU_DIR

# launcher enter to scratch dir with copied:
# geom.in, input.in ABINITIO folder etc. which are read by init.py
cwd = os.getcwd()
sys.path.append(cwd)

modules_files = [
    'init.py',
    'masses.py',
    'errors.py',
    'prints.py',
    'propagates.py',
]

try:
    SNAFU_EXE = os.environ['SNAFU_DIR']
    sys.path.append(os.path.join(SNAFU_EXE,"snafu"))
    sys.path.append(SNAFU_EXE)
    from snafu.init import file_check, read_input
    from snafu.init import read_geoms, read_velocs
    from snafu.init import create_output_file, init_forces_potenergs
    from snafu.init import com_removal
    from snafu.masses import assign_masses
    from snafu.errors import error_exit
    from snafu.prints import print_positions, print_velocities
    from snafu.prints import print_snafu
    from snafu.propagates import update_velocities, update_positions
    from snafu.propagates import calc_forces, calc_energies
except ImportError as ime: # some module file is missing or was renamed
    print("Module {} not found.".format(ime.name),
          "Make sure that {} dir contains snafu folder".format(SNAFU_EXE),
          "contains: {}".format('\n'.join(modules_files)))
    exit(1)
except KeyError as ke: # if SNADU_DIR is not exported into env vars
    print("SNAFU_DIR is not set.",
          "See '.bashrc' file in your home directory", 
          "or use 'env' command and make sure $SNAFU_DIR is exported.",
          "\nHint: export SNAFU_DIR=/path/to/SNAFU")
    exit(1)
else:
    print("All modules loaded succesfully. Starting...\n \n")


# CONVERSION CONSTANTS
au_fs = 0.02418884326505e0      #atomic units to femtosecs
au_eV = 27.21139
amu   = 1822.888484264545e0     # atomic mass unit  me = 1 AMU*atomic weight
ang_bohr = 1.889726132873e0     # agstroms to bohrs

# RUN VARIABLES
step = 0
dE = 0.0   # energy change since initial energies
Etot_init = 0.0  # setting variable , total energy at the beginning 
time = 0.0
liner = ("_")*70
debug = 1   
#----------------INIT-------------------------------------------------------------------------------

if __name__ == "__main__":
    print_snafu()
    startTime = datetime.now()
    
    print("Simulation started at: {}\n".format(startTime),
          "Python base: {}".format(sys.base_exec_prefix),
          "version: {}\n".format(sys.version[:5]),
          "System platform: {}\n".format(sys.platform),
          "Running executable: {}\n".format(sys.path[-1])
          )
    # local runs doent creat HOST env var, qsub SGE system does
    try: 
        print("Working directory: {} on {}.".format(cwd,os.environ['HOSTNAME']))
    except KeyError: 
        print("Working directory: {}".format(cwd))
    print(liner)
    
# FILE CHECK - OBTAIN PATH TO FILES
    input_file_path,geom_file_path,veloc_file_path,veloc_init = file_check(cwd)
    
#READ INPUT VARIABLE (read as strings) - SET THEM AS GLOBAL VARIABLES: 
 
    input_vars, ab_initio_file_path = read_input(cwd,input_file_path)
    globals().update(input_vars)  
    natoms   = int(natoms)                                    
    maxsteps = int(maxsteps)
    state = int(init_state)   # initial or restart state
    dt = float(timestep)      
    nstates = int(nstates)
    
#READ INITIAL GEOMETRY AND VELOCITIES AND CREATE ARRAYS FOR FORCES:
    
    # names and positions array (1D)
    at_names, x, y, z = read_geoms(natoms,geom_file_path)    
                           
    # Velocity array (1D) 
    vx, vy, vz = read_velocs(veloc_init,natoms,veloc_file_path)                    
    
    # Empty forces and potential energy arrays (1D)
    fx, fy, fz, fx_new, fy_new,fz_new, pot_eners = init_forces_potenergs(
        natoms,nstates)         
    
#OBTAIN MASSES:
    print(liner)
    masses = assign_masses(at_names) 
    am = [ mm *  amu for mm in masses]  # atomic mass units conversion    
    print("At  X     Y     Z   MASS:")
    for iat in range(0,natoms):
        print("".join("%2s" " " "%3.3f"  %(at_names[iat], x[iat])),
              " %3.6f %2.6f %2.6f"  %(y[iat],z[iat],masses[iat]))  
    print(liner)     
    
#CREATE OUTPUT FILES: 
   # TO DO file check, if these files exist and not restart then error!!!!
   # where to store propagated position, velocities, observables
    # geom.dat hold current geometry for which to compute E, grads  
    #files = [ "energies.dat", "velocities.xyz", "gradients.dat", "movie.xyz" ]  
    #create_output_file(files)

#---------------INIT DONE-------------------------------------------------------------------    

# CENTER OF MASS REMOVAL
    x, y, z = com_removal(x,y,z,am)
    
# CALC INITIAL ENERGIES AND GRADIENTS
    print("Step      Time/fs     Energy change from start/eV  Hoppping") 
    
    fx, fy, fz, pot_eners = calc_forces(step, natoms, at_names, state,
        nstates, ab_initio_file_path,x, y, z, fx, fy, fz, pot_eners)  
    
    Ekin, Epot, Etot, dE = calc_energies(step, time, natoms, am, state, 
        pot_eners, vx, vy, vz, Etot_init) #Etot_init = 0
    
    # Total energy at the beginning to calc. energy changes during propagation    
    Etot_init = Etot  
    
    
# MAIN LOOP 
    for step in range(1,maxsteps+1):
        # new positions (t+dt)
        x, y, z = update_positions(natoms,dt,am,
            x,y,z,vx,vy,vz,fx,fy,fz) 
        
        # calc forces for new positions       
        fx_new, fy_new, fz_new, pot_eners = calc_forces(
            step,natoms, at_names, state, nstates,ab_initio_file_path,
            x, y, z, fx_new, fy_new, fz_new,pot_eners)  
       
        # propagate velocities using new forces
        vx, vy, vz = update_velocities(natoms,dt,am,vx,vy,vz,fx,fy,fz,
            fx_new,fy_new,fz_new) 
        
        # copied list instead of just referencing | or slice it
        fx = np.copy(fx_new) 
        fy = np.copy(fy_new)
        fz = np.copy(fz_new)
        
        #if hopping == "1":
          #alcc_hop
        # vel_adjustment
        
        hop = "No"
        time = step * dt*  au_fs 
        Ekin, Epot, Etot, dE = calc_energies(step, time, natoms, 
            am, state, pot_eners, vx, vy, vz, Etot_init)
        
        print(" {:<3d} {:>10.2f} {:>20.4e}".format(step,time,dE * au_eV),
              " {:>20s}".format(hop))
        
        # save positions and velocities
        print_positions(step,time,natoms, at_names, x, y, z)   
        print_velocities(step,time,natoms, at_names, vx, vy, vz)
      
    
    print("JOB completed.") 
    print(liner)
    
    stopTime = datetime.now()
    simtime = (datetime.now() - startTime)
    print("Simulation ended at: {}".format(stopTime))
    print("Overall simulation time (hh:mm:ss): {}".format(simtime))     

