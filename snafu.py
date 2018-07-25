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
np.set_printoptions(precision=8)  # for print and stability testing

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
]

# if some of the imports dont work in other modules
# ImportError will be raised with None name
try:
    SNAFU_EXE = os.environ['SNAFU_DIR']
    sys.path.append(os.path.join(SNAFU_EXE, "snafu"))
    sys.path.append(SNAFU_EXE)
    from inits import (
        file_check, read_input, read_geoms, read_velocs,
        create_output_file, com_removal, init_forces_energs
    )
    from masses import assign_masses
    from errors import error_exit
    from prints import (
        print_positions, print_velocities,
        print_snafu, print_energies
    )
    from propagates import (
        calc_forces, calc_energies,
        update_velocities, update_positions
    )
except ImportError as ime:
    if ime.name is None:
        print("Import in some of the modules in snadu dir failed")
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

# CONVERSION CONSTANTS
au_fs = 0.02418884326505e0   # atomic units to femtosecs
au_eV = 27.21139
amu = 1822.888484264545e0    # atomic mass unit me = 1 AMU*atomic weight
ang_bohr = 1.889726132873e0  # agstroms to bohrs

# RUN VARIABLES
step = 0
dE = 0.0  # energy change since initial energies
Etot_init = 0.0  # setting variable , total energy at the beginning
time = 0.0
liner = ("_") * 70

# ---------------INIT---------------------------------------------------------

if __name__ == "__main__":

    print("Simulation started at: {}\n".format(startTime),
          "Python base: {}".format(sys.base_exec_prefix),
          "version: {}\n".format(sys.version[:5]),
          "System platform: {}\n".format(sys.platform),
          "Running executable: {}\n".format(sys.path[-1])
          )
    # local runs dont create HOST env var, qsub SGE system does
    try:
        print("Working directory: {} on {}.".format(cwd,
                                                    os.environ['HOSTNAME']))
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

    at_names, x, y, z = read_geoms(natoms, geom_file_path)

    vx, vy, vz = read_velocs(init_vel, natoms, vel_file_path)

    fx, fy, fz, fx_new, fy_new, fz_new, pot_eners = init_forces_energs(
        natoms, nstates)

    # OBTAIN ATOMIC MASSES:
    masses = assign_masses(at_names)
    am = [mm * amu for mm in masses]  # atomic mass units conversion
    print("At  X     Y     Z   MASS:")
    for iat in range(0, natoms):
        print("".join("%2s" " " "%3.3f" % (at_names[iat], x[iat])),
              " %3.6f %2.6f %2.6f" % (y[iat], z[iat], masses[iat]))
    print(liner)

    # TO DO CHECK RESTART files OR if files exist and not restart then error!!

    # CENTER OF MASS REMOVAL
    x, y, z = com_removal(x, y, z, am)

    # CALC INITIAL ENERGIES AND GRADIENTS
    print("Step      Time/fs     Energy change from start/eV  Hoppping")

    fx, fy, fz, pot_eners = calc_forces(step, at_names, state, nstates,
                                        x, y, z, fx, fy, fz, pot_eners,
                                        ab_initio_file_path)

    Ekin, Epot, Etot, dE = calc_energies(step, time, natoms, am, state,
                                         pot_eners, vx, vy, vz, Etot_init)
    Etot_init = Etot

    # MAIN LOOP
    for step in range(1, maxsteps + 1):

        x, y, z = update_positions(dt, am, x, y, z, vx, vy, vz, fx, fy, fz)

        fx_new, fy_new, fz_new, pot_eners = calc_forces(step, at_names, state,
                                                        nstates, x, y, z,
                                                        fx_new, fy_new, fz_new,
                                                        pot_eners,
                                                        ab_initio_file_path)

        vx, vy, vz = update_velocities(dt, am, vx, vy, vz, fx, fy, fz,
                                       fx_new, fy_new, fz_new)

        # copied list instead of just referencing | or slice it
        fx = np.copy(fx_new)
        fy = np.copy(fy_new)
        fz = np.copy(fz_new)

        # if hopping == "1": vel_adjustment

        hop = "No"
        time = step * dt * au_fs
        Ekin, Epot, Etot, dE = calc_energies(step, time, natoms, am,
                                             state, pot_eners,
                                             vx, vy, vz, Etot_init)

        print(" {:<3d} {:>10.2f} {:>20.4e}".format(step, time, dE * au_eV),
              " {:>20s}".format(hop))

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
