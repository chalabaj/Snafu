# SYSTEM IMPORTS:
import math
import sys
import os
import random
import time
import re
import numpy as np
from datetime import datetime
startTime = datetime.now()


# ENVIRONMENT LAYER & LOCAL IMPORT
# find local modules
# SNAFU_DIR=/path/to/SNAFU/snafu variable must be set before launch
# either export directly or in .bashrc SNAFU_DIR
# launcher enter to scratch dir with copied input files
cwd = os.getcwd()
sys.path.append(cwd)

modules_files = [
    'inits.py',
    'masses.py',
    'errors.py',
    'prints.py',
    'propagates.py',
    'landauzener.py',
    'constants.py'
]

# if some of the imports fail in modules ImportError rais None name
try:
    SNAFU_EXE = os.environ['SNAFU_DIR']
    sys.path.append(os.path.join(SNAFU_EXE, "snafu"))
    sys.path.append(SNAFU_EXE)
    from inits import (
        file_check, read_input, read_geoms, read_velocs,
        com_removal, init_forces, init_energies
    )
    from masses import assign_masses
    from errors import error_exit
    from prints import (
        print_positions, print_velocities, print_snafu
    )
    from propagates import (
        calc_forces, calc_energies,
        update_velocities, update_positions, rescale_velocities,
        adjust_velocities
    )
    from landauzener import (
        calc_hopp
    )
    from constants import *
except ImportError as ime:
    # module could have been removed or module file renamed
    if ime.name is None:  
        print("Import in some of the modules ({})".format(ime.name),
              "in snafu dir failed. Exiting...")
        exit(1)
    else:
        print("Module {} not found.".format(ime.name),
              "Make sure that {} contains snafu folder".format(SNAFU_EXE),
              "with: {}".format('\n'.join(modules_files)))
        exit(1)
except KeyError as ke:
    print("SNAFU_DIR is not set.",
          "See '.bashrc' file in your home directory",
          "or use 'env' command and make sure $SNAFU_DIR is exported.",
          "\nHint: export SNAFU_DIR=/path/to/SNAFU")
    exit(1)
else:
    print("All modules loaded succesfully. Starting...\n \n")
    print_snafu()


liner = ("_") * 100

# ---------------INIT---------------------------------------------------------

if __name__ == "__main__":

    print("Simulation started at: {}".format(startTime),
          "\nPython base: {}".format(sys.base_exec_prefix),
          "version: {}".format(sys.version[:5]),
          "\nSystem platform: {}".format(sys.platform),
          "\nRunning executable: {}".format(sys.path[-1])
          )
    # local runs dont create HOST env var, qsub SGE system does
    try:
        print("Working directory: {}".format(cwd),
              "on {}.".format(os.environ['HOSTNAME'])
             )
    except KeyError:
        print("Working directory: {}".format(cwd))
    print(liner)

    # FILE CHECK - OBTAIN PATH TO FILES
    input_file_path, geom_file_path, vel_file_path, init_vel = file_check(cwd)

    # READ INPUT VARIABLES SET THEM AS GLOBAL VARIABLES
    input_vars, ab_initio_file_path = read_input(cwd, input_file_path)
    globals().update(input_vars)
    natoms = int(natoms)
    maxsteps = int(maxsteps)
    state = int(init_state)   # initial or restart state
    dt = float(timestep)
    nstates = int(nstates)
    print(liner)

    # READ INITIAL GEOMETRY AND VELOCITIES AND CREATE ARRAYS FOR FORCES
    # (1D array)
    hop = False
    step = 0
    
    at_names, x, y, z, x_new, y_new, z_new = read_geoms(natoms,
                                                        geom_file_path)

    vx, vy, vz = read_velocs(init_vel, natoms, vel_file_path)

    fx, fy, fz, fx_new, fy_new, fz_new = init_forces(natoms, nstates)

    pot_eners = init_energies(nstates)
    
    # OBTAIN ATOMIC MASSES:
    masses = assign_masses(at_names)
    am = [mm * AMU for mm in masses]  # atomic mass units conversion
    print("At    X       Y         Z         MASS:")
    for iat in range(0, natoms):
        print("".join("%2s" " " "%3.3f" % (at_names[iat], x[iat])),
              " %3.6f %2.6f %2.6f" % (y[iat], z[iat], masses[iat]))
    print(liner)

    # TO DO CHECK RESTART files OR if files exist and not restart then error!!

    # CENTER OF MASS REMOVAL
    x, y, z = com_removal(x, y, z, am)

    # CALC INITIAL ENERGIES AND GRADIENTS

    fx, fy, fz, pot_eners = calc_forces(step, at_names, state, nstates,
                                        x, y, z, fx, fy, fz, pot_eners,
                                        ab_initio_file_path)
    pot_eners_array = np.copy(pot_eners)                                    
    dE = 0.0  # energy change since initial energies
    Etot_init = 0.0  # setting variable , total energy at the beginning
    Etot_prev = 0.0
    time = 0.0
    prob = 0.0
    Ekin, Epot, Etot, dE, dE_step = calc_energies(step, time, natoms, am,
                                                  state, pot_eners,
                                                  vx, vy, vz, Etot_init,
                                                  Etot_prev)
    Etot_init = Etot
    
    print("Step    Time/fs      dE_drift/eV      dE_step/eV      Hoppping",
          "  State  Max. prob")
    # MAIN LOOP
    for step in range(1, maxsteps + 1):
          
        x_new, y_new, z_new = update_positions(dt, am, 
                                               x, y, z, 
                                               vx, vy, vz, 
                                               fx, fy, fz)

        fx_new, fy_new, fz_new, pot_eners = calc_forces(step, at_names, 
                                                        state, nstates, 
                                                        x_new, y_new, z_new,
                                                        fx_new, fy_new,fz_new,
                                                        pot_eners,
                                                        ab_initio_file_path)
        
        if not method == "bomd":
            print("evaluate hop")
            if step >= 2:
                hop, outstate, v_scal_fac, prob = calc_hopp(method, state, 
                                                        pot_eners, 
                                                        pot_eners_array, 
                                                        Ekin, dt)

                if hop:
                    state = outstate
               
                    # use R  from prev. step, hop, cacl f new state
                    fx_new, fy_new, fz_new, pot_eners = calc_forces(
                        step, at_names, state, nstates, 
                        x, y, z,
                        fx_new, fy_new, fz_new, 
                        pot_eners, ab_initio_file_path)
                   #simple scaling or updatre velocities with new state forces
                    if not int(vel_adj):
                        print("Rescale")
                        vx, vy, vz = rescale_velocities(vx, vy, vz, v_scal_fac)
                    else:
                        vx, vy, vz = adjust_velocities(dt, am,
                                                      vx, vy, vz,
                                                      fx, fy, fz,
                                                      fx_new, fy_new, fz_new) 
                    EE = 0
                    for iat in range(0,natoms):
                        vvv = vx[iat] ** 2 + vy[iat] ** 2 + vz[iat] ** 2
                        EE = EE + (0.5 * am[iat] * vvv)
                    print("Scaled Ekin: {}, old Ekin {}".format(EE, Ekin))

                    # now finish the propagation step on new PES
                    x_new, y_new, z_new = update_positions(dt, am, 
                                                           x, y, z, 
                                                           vx, vy, vz, 
                                                           fx_new, fy_new,
                                                           fz_new)

                    fx_new, fy_new, fz_new, pot_eners = calc_forces(
                        step, at_names, state, nstates, 
                        x_new, y_new, z_new,
                        fx_new, fy_new, fz_new,
                        pot_eners, ab_initio_file_path)

                pot_eners_array = np.delete(pot_eners_array, 0, axis = 0)
                pot_eners_array = np.vstack((pot_eners_array, pot_eners)) 

            else:
                pot_eners_array = np.vstack((pot_eners_array, pot_eners)) 
            
        vx, vy, vz = update_velocities(dt, am, 
                                       vx, vy, vz,
                                       fx, fy, fz,
                                       fx_new, fy_new, fz_new)
        fx = np.copy(fx_new)
        fy = np.copy(fy_new)
        fz = np.copy(fz_new)
        x = np.copy(x_new)
        y = np.copy(y_new)
        z = np.copy(z_new)

        time = step * dt * AU_FS
        Etot_prev = Etot
        Ekin, Epot, Etot, dE, dE_step = calc_energies(step, time, natoms, am,
                                                      state, pot_eners,
                                                      vx, vy, vz, Etot_init,
                                                      Etot_prev)
        if hop:
            print("new X Ekin {}".format(Ekin))
        print(" {:<4d} {:<12.2f} {:>12.4e}".format(step, time, dE * AU_EV),
              " {:>12.4f}".format(dE_step*AU_EV),
              " {:>12s}     {}        {} \n".format(str(hop), state, prob))

        # save positions and velocities
        print_positions(step, time, natoms, at_names, x, y, z)
        print_velocities(step, time, natoms, at_names, vx, vy, vz)

    print("JOB completed.")
    print(liner)
    stopTime = datetime.now()
    simtime = (datetime.now() - startTime)
    print("Simulation ended at: {}".format(stopTime))
    print("Overall simulation time (hh:mm:ss): {}".format(simtime))
exit(0)
