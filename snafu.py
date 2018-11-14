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
    'constants.py',
    'restart.py'
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
        print_positions, print_velocities, print_snafu,
        print_state
    )
    from propagates import (
        calc_forces, calc_energies,
        update_velocities, update_positions, rescale_velocities,
        adjust_velocities
    )
    from landauzener import (
        calc_hopp
    )
    from restarts import (
        print_restart, check_restart, read_restart
    )
    from constants import *
    from tera_propagates import (
        finish_tera, exit_tera, tera_connect, tera_init

    )
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

# ---------------INIT--------------------------------------------------------

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
    # DEFAULTS: 
    # TODO: move these to init routine and import them before reading the input
    ener_thresh = 1.000
    hop_thresh = 0.5
    dt = 4.00
    hop = False
    step = 0
    dE = 0.0  # energy change since initial energies
    Etot_init = 0.0  # setting variable , total energy at the beginning
    Etot_prev = 0.0
    sim_time = 0.0
    prob = 0.0
    tera_mpi = 0

    input_vars, ab_initio_file_path = read_input(cwd, input_file_path)
    print(liner)
    globals().update(input_vars)
    try:
        natoms = int(natoms)
        maxsteps = int(maxsteps)
        state = int(init_state)   # initial or restart state
        dt = float(timestep)
        nstates = int(nstates)
        ener_thresh = float(ener_thresh)
        hop_thresh = float(hop_thresh)
        vel_adj = int(vel_adj)
        restart = int(restart)
        restart_write = int(restart_write)
        tera_mpi = int(tera_mpi)
    except ValueError as VE:
        print(VE)
        error_exit(9, " ")
    
    pot_eners = init_energies(nstates)
    fx, fy, fz, fx_new, fy_new, fz_new = init_forces(natoms, nstates)
    # READ INITIAL GEOMETRY AND VELOCITIES AND CREATE ARRAYS FOR FORCES
    rst_file_path = check_restart(restart, cwd)
    print(liner)
    
    if tera_mpi:
        byte_coords = np.dstack((x,y,z)))  # x,y,z must stack to single array
        comm = tera_connect()
        MO, MO_old, CiVecs, CiVecs_old, NAC, blob, \
        blob_old, SMatrix, civec_size, nbf_size,  \
        blob_size, TDip, Dip = tera_init(comm, at_names, natoms, nstates,
                                         byte_coords) 
    print(liner)

    if restart == 0:
    
        at_names, x, y, z, x_new, y_new, z_new = read_geoms(natoms,
                                                            geom_file_path)
        vx, vy, vz = read_velocs(init_vel, natoms, vel_file_path)

        # OBTAIN ATOMIC MASSES:
        masses = assign_masses(at_names)
        am = [mm * AMU for mm in masses]  # atomic mass units conversion
        # CENTER OF MASS REMOVAL 
        x, y, z = com_removal(x, y, z, am)
        # CALC INITIAL ENERGIES AND GRADIENTS
        fx, fy, fz, pot_eners = calc_forces(step, at_names, state, nstates,
                                            x, y, z, fx, fy, fz, pot_eners,
                                            ab_initio_file_path)
        pot_eners_array = np.copy(pot_eners)      
        
        Ekin, Epot, Etot, dE, dE_step = calc_energies(step, sim_time, natoms, am,
                                                      state, pot_eners,
                                                      vx, vy, vz, Etot_init,
                                                      Etot_prev, ener_thresh,
                                                      restart)
        Etot_init = Etot
        init_step = 1
    else:
        init_step, at_names, state, \
        x, y, z, vx, vy, vz, fx, fy, fz, \ 
        Ekin, Epot, Etot, Etot_init, \
        pot_eners_array = read_restart(rst_file_path, natoms)
        
        masses = assign_masses(at_names)
        am = [mm * AMU for mm in masses]  # atomic mass units conversion
        x_new = np.zeros_like(x)
        y_new = np.zeros_like(y)
        z_new = np.zeros_like(z)
        fx_new = np.zeros_like(fx)
        fy_new = np.zeros_like(fy)
        fz_new = np.zeros_like(fz)
        init_step = init_step + 1
    print("Initial geometry:\n",
          "At    X         Y          Z         MASS:")
    xx = [xxx*BOHR_ANG for xxx in x.tolist()]
    yy = [yyy*BOHR_ANG for yyy in y.tolist()]
    zz = [zzz*BOHR_ANG for zzz in z.tolist()]
    for iat in range(0, natoms):
        print("{} {:12.8f}".format(at_names[iat], xx[iat]),
              "{:12.8f} {:12.8f}".format(yy[iat], zz[iat]),
              "{:12.8f}".format(masses[iat]))

    print("Initial velocities:\n",
          "At    VX       VY       VZ         ")
    for iat in range(0, natoms):
        print("".join("%2s" " " "%3.3f" % (at_names[iat], vx[iat])),
              " %3.6f %2.6f " % (vy[iat], vz[iat]))
    print("{}".format(liner),
          "\nStep    Time/fs  dE_drift/eV   dE_step/eV    Hop  State") 

    #-------------------MAIN LOOP-----------------------------------------
    for step in range(init_step, maxsteps + 1):

        x_new, y_new, z_new = update_positions(dt, am, 
                                               x, y, z,
                                               x_new, y_new, z_new, 
                                               vx, vy, vz, 
                                               fx, fy, fz)

        fx_new, fy_new, fz_new, pot_eners = calc_forces(step, at_names, 
                                                        state, nstates, 
                                                        x_new, y_new, z_new,
                                                        fx_new, fy_new, fz_new,
                                                        pot_eners,
                                                        ab_initio_file_path)

        if not method == "bomd":
            if step >= 2:
                hop, outstate, v_scal_fac, prob = calc_hopp(method, state,
                                                            pot_eners, 
                                                            pot_eners_array,
                                                            Ekin, dt,
                                                            hop_thresh)

                if hop:
                    state = outstate
                    # use XYZ from prev. step to cacl F for a new state
                    
                    fx_new, fy_new, fz_new, pot_eners = calc_forces(
                        step, at_names, state, nstates, 
                        x, y, z,
                        fx_new, fy_new, fz_new, 
                        pot_eners, ab_initio_file_path)

                    #simple scaling or updatre velocities with new state forces

                    if not vel_adj:
                        vx, vy, vz = rescale_velocities(vx, vy, vz, v_scal_fac)
                    else:
                        vx, vy, vz = adjust_velocities(dt, am,
                                                       vx, vy, vz,
                                                       fx, fy, fz,
                                                       fx_new, fy_new, fz_new) 

                    # FXFYFZ for  XYZ(t = hop) in the new state already
                    fx = np.copy(fx_new)
                    fy = np.copy(fy_new)
                    fz = np.copy(fz_new)
                    EE = 0
                    for iat in range(0,natoms):
                        vvv = vx[iat] ** 2 + vy[iat] ** 2 + vz[iat] ** 2
                        EE = EE + (0.5 * am[iat] * vvv)
                    print("Old Ekin: {} \nScaled/Adjusted Ekin {}\n.".format(Ekin, EE))

                    # now finish the propagation step on new PES

                    x_new, y_new, z_new = update_positions(dt, am, 
                                                           x, y, z, 
                                                           x_new, y_new, z_new, 
                                                           vx, vy, vz, 
                                                           fx_new, fy_new,
                                                           fz_new)


                    fx_new, fy_new, fz_new, pot_eners = calc_forces(
                        step, at_names, state, nstates, 
                        x_new, y_new, z_new, fx_new, fy_new, fz_new,
                        pot_eners, ab_initio_file_path)

                pot_eners_array = np.delete(pot_eners_array, 0, axis = 0)
                pot_eners_array = np.vstack((pot_eners_array, pot_eners))  #  keep last two steps

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

        sim_time = step * dt * AU_FS
        Etot_prev = Etot
        Ekin, Epot, Etot, dE, dE_step = calc_energies(step, sim_time, natoms, am,
                                                      state, pot_eners,
                                                      vx, vy, vz, Etot_init,
                                                      Etot_prev, ener_thresh,
                                                      restart)

        # print("Ekin {}, Epot {}, Etot{}".format(Ekin, Epot, Etot))
        print(" {:<6d}  {:<7.4f}  {:<12.4f}".format(step, sim_time, dE * AU_EV),
              " {:<13.4f}".format(dE_step * AU_EV),
              "{}     {}".format(str(hop)[0], state))
        #print("-----------------------------------------------------")

        # SAVE POSITION AND VELOCITIES AND RESTART
        print_positions(step, sim_time, natoms, at_names, x, y, z, restart)
        print_velocities(step, sim_time, natoms, at_names, vx, vy, vz, restart)
        print_state(step, sim_time, state, restart)

        print_restart(step, sim_time, natoms, at_names, state, timestep,
                      x, y, z, vx, vy, vz, fx, fy, fz,
                      Ekin, Epot, Etot, Etot_init, pot_eners_array,
                      restart_write)
    # FINAL PRINTS
    print(liner)
    print("#####JOB DONE.############")
    print("See output files:",
          "\nmovie.xyz, velocities.xyz,\nPES.dat, energies.dat,\nstate.dat")
    print(liner)
    stopTime = datetime.now()
    calc_time = (datetime.now() - startTime)
    print("Simulation ended at: {}".format(stopTime))
    print("Overall simulation sim_time (hh:mm:ss): {}".format(calc_time))
exit(0)
