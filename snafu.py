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

modules_files = ['masses.py',
'constants.py',
'landauzener.py',
'restarts.py',
'propagates.py',
'prints.py',
'errors.py',
'defaults.py',
'inits.py',
'tera-propagate.py'
]

# if some of the imports fail during in-modules import ImportError rais None name
try:
    SNAFU_EXE = os.environ['SNAFU_DIR']
    sys.path.append(os.path.join(SNAFU_EXE, "snafu"))
    sys.path.append(SNAFU_EXE)
    from inits import (
        file_check, read_input, check_output_file, read_geoms, read_velocs,
        com_removal, init_fep_arrays
    )
    from masses import assign_masses
    from errors import error_exit
    from prints import (
        print_positions, print_velocities, print_snafu,
        print_state, print_energies, print_pes
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
        print_restart, check_restart_files, read_restart
    )
    from constants import *   #  import conversion factors and default values
    
    from defaults import *    # import all defualt values, only here otherwise could overwritte in some modules
    #from tera_propagates import (
    #    finish_tera, exit_tera, tera_connect, tera_init
    #)
except ImportError as ime:
    # module could have been removed or different module name, e.g. renamed in module file
    if ime.name is None:  
        print("Import in some modules {}".format(ime),
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
    print("All modules loaded succesfully.\n \n")
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
              "on {}.".format(os.environ['HOSTNAME']))
    except KeyError:
        print("Working directory: {}".format(cwd))
    print(liner)

    # FILE CHECK - OBTAIN PATH TO FILES
    input_file_path, geom_file_path, vel_file_path, init_vel = file_check(cwd)
    print(liner)

    #  READ INPUT OPTIONS AND SET THEM AS VARIABLES:
    input_vars, ab_initio_file_path = read_input(cwd, input_file_path)
    #  Need to pass variables to function so that modules wont used defaults from constats mod
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
        restart_freq = int(restart_freq)
        tera_mpi = int(tera_mpi)
        write_freq = int(write_freq)            
        print("Simulation will statr with the following parameters:\n",
              "{} = {}\n".format("natoms", natoms),
              "{} = {}\n".format("maxsteps",maxsteps),
              "{} = {}\n".format("initial state", state),
              "{} = {}\n".format("timeste", dt),
              "{} = {}\n".format("nstates",nstates),
              "{} = {}\n".format("ener_thresh", ener_thresh),
              "{} = {}\n".format("hop_thresh", hop_thresh),
              "{} = {}\n".format("vel_adj", vel_adj),
              "{} = {}\n".format("restart", restart),
              "{} = {}\n".format("restart_freq", restart_freq),
              "{} = {}\n".format("tera_mpi", tera_mpi),
              "{} = {}\n".format("write_freq",write_freq),
              "{} = {}\n".format("abinitio file",ab_initio_file_path),
              "{} = {}\n".format("method", method))
    except ValueError as VE:
        error_exit(9, str(VE))

    fx, fy, fz, fx_new, fy_new, fz_new, \
    pot_eners, x_new, y_new, z_new = init_fep_arrays(natoms, nstates)

    # READ INITIAL OR RESTART DATA
    rst_file_path = check_restart_files(restart, cwd)
    print(liner)
    if restart == 0:
        at_names, x, y, z  = read_geoms(natoms, geom_file_path)
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
                                                      Etot_prev, ener_thresh)
        Etot_init = Etot
        init_step = 1
    else:
        init_step, at_names, state, \
        x, y, z, vx, vy, vz, fx, fy, fz, \
        Ekin, Epot, Etot, Etot_init, \
        pot_eners_array = read_restart(rst_file_path, natoms)
        
        masses = assign_masses(at_names)
        am = [mm * AMU for mm in masses]  # atomic mass units conversion
        init_step = init_step + 1

    check_output_file(cwd, natoms, restart, init_step)
    
    print("Initial geometry:\n",
          "At    X         Y         Z         MASS:")
    xx = (x*BOHR_ANG).tolist()
    yy = (y*BOHR_ANG).tolist()
    zz = (z*BOHR_ANG).tolist()
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
    
    #if tera_mpi:
    #    byte_coords = np.dstack((x,y,z))  # x,y,z must stack to single array
    #    comm = tera_connect()
    #    MO, MO_old, CiVecs, CiVecs_old, NAC, blob, \
    #    blob_old, SMatrix, civec_size, nbf_size,  \
    #    blob_size, TDip, Dip = tera_init(comm, at_names, natoms, nstates,
    #                                     byte_coords) 
    print(liner)
    with open('movie.xyz', 'a') as mov_file, \
         open('energies.dat', 'a') as eners_file, \
         open('PES.dat', 'a') as pes_file, \
         open('velocities.dat', 'a') as vel_file, \
         open('state.dat', 'a') as state_file, \
         open('restart.in', 'w') as rsf_file:

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
                                                          Etot_prev, ener_thresh)
    
            # print("Ekin {}, Epot {}, Etot{}".format(Ekin, Epot, Etot))
            print(" {:<6d}  {:<7.4f}  {:<12.4f}".format(step, sim_time, dE * AU_EV),
                  " {:<13.4f}".format(dE_step * AU_EV),
                  "{}     {}".format(str(hop)[0], state))
            # SAVE POSITION, VELOCITIES, ENERGIES AND RESTART
            if (step%write_freq == 0):
                print_energies(step,write_freq, sim_time, Ekin, Epot, Etot, dE, dE_step, eners_file)
                print_pes(step, write_freq, sim_time, pot_eners, pes_file)
                print_positions(step, sim_time, natoms, at_names, x, y, z, mov_file)
                print_velocities(step, sim_time, natoms, at_names, vx, vy, vz, vel_file)
                print_state(step, write_freq, sim_time, state, state_file)
            
            print_restart(step, sim_time, natoms, at_names, state, timestep,
                          x, y, z, vx, vy, vz, fx, fy, fz,
                          Ekin, Epot, Etot, Etot_init, pot_eners_array, restart_freq, rsf_file)
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
