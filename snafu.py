"""
This is main file which runs the LZ simulations.
Import required module.
Read inputs (input option, geometry, velocities, restarts)
Run Velocity verlet integrator MD with LZ hopping algorithm at each step
"""

# SYSTEM IMPORTS:
import sys
import os
import numpy as np
from datetime import datetime
startTime = datetime.now()

# ENVIRONMENT LAYER 
try:
    cwd = os.getcwd()
    sys.path.append(cwd)
    SNAFU_EXE = os.environ['SNAFU_DIR']
    sys.path.append(os.path.join(SNAFU_EXE, "snafu"))
    sys.path.append(SNAFU_EXE)
except KeyError as ke:
    # SNAFU_DIR=/path/to/SNAFU/snafu variable must be set before launch
    print("SNAFU_DIR is not set.",
          "See '.bashrc' file in your home directory",
          "or use 'env' command and make sure $SNAFU_DIR is exported.",
          "\nHint: export SNAFU_DIR=/path/to/SNAFU")
    sys.exit(1)

# TERA MPI INTERFACE SMOOTH EXIT:    
# load excepthook as soon as possible in order to prevent ancaught exceptions which can cause MPI deadlock state
# This will mainly check code and runtime errors, wrong inputs; known errors are handled in error_exit module
# sys.excepthook has some collisions and does NOT work if syntax error is in some imported module bellow, use mangling
try: 
    tera_mpi = int(os.environ['MPI_TERA'])
    if tera_mpi:
        from tera_propagates import (finish_tera, tera_connect, tera_init, global_except_hook)
        sys.excepthook = global_except_hook    
except KeyError as ke:
     print("MPI_TERA variable was not exported, assuming MPI_TERA=0. Warning: this may cause deadlock if MPI has been already initiated")
     tera_mpi = 0

# LOCAL IMPORT OF SNAFU MODULES 
try:
    modules_files = ['masses.py','constants.py','landauzener.py','restarts.py','propagates.py','prints.py',
                     'errors.py','defaults.py','inits.py','tera_propagates.py']    
    from inits import (
        file_check, read_input, 
        check_output_file, read_geoms, 
        read_velocs, com_removal, 
        init_fep_arrays
    )
    from masses import assign_masses
    from errors import error_exit
    from prints import (
        print_positions, print_velocities, 
        print_snafu, print_state, 
        print_energies, print_pes
    )
    from propagates import (
        calc_forces, calc_energies,
        update_velocities, update_positions, 
        rescale_velocities, adjust_velocities
    )
    from landauzener import (
        calc_hopp
    )
    from restarts import (
        print_restart, check_restart_files, read_restart
    )
    from constants import *   #  conversion factors (e.g. BOHR to ANG units)
    from defaults import *    #  All default values for the code, import only here
except ImportError as ime:
    # if some of the following imports fail during in-modules-import: ImportError raise None
    # module could have been removed or different module name, e.g. renamed in module file
    if ime.name is None:  
        error_exit(19, "Import or internal error in some module {}.".format(ime))
    else:
        print("Module {} not found!!!".format(ime.name),
              "Make sure that {} contains snafu folder".format(SNAFU_EXE),
              "with: {}.".format('\n'.join(modules_files)),
              "\nOr check import in the wrong module")
        error_exit(19, "Some python file probably missing {}".format(ime))
else:
    print_snafu()
    print("\nAll modules loaded.")

# ---------------INITIALIZATION PART--------------------------------------------------------
if __name__ == "__main__":
    print("Simulation started at: {}".format(startTime),
          "\nPython: {}".format(sys.base_exec_prefix),
          "version: {}".format(sys.version[:5]),
          # "\nSystem platform: {}".format(sys.platform),
          "\nSystem path: {}".format(sys.path[0]),
          "\n".join(sys.path[1:]) 
          )
    # Local runs dont create HOSTNAME variable, qsub SGE system does
    try:
        print("Working directory: {}".format(cwd),
              "on {}.".format(os.environ['HOSTNAME']))
    except KeyError:
        print("Working directory: {}".format(cwd))
    print(liner)

    # FILE CHECK - OBTAIN PATHS TO ALL FILES
    input_file_path, geom_file_path, vel_file_path, init_vel = file_check(cwd)

    #  READ INPUT OPTIONS AND SET THEM AS VARIABLES:
    # with exception of natoms, nstates, following option can be changed by restart
    input_vars, ab_initio_file_path = read_input(cwd, input_file_path)
    globals().update(input_vars)
    try:
        natoms = int(natoms)
        maxsteps = int(maxsteps)
        state = int(init_state)   
        dt = float(timestep)
        nstates = int(nstates)
        ener_thresh = float(ener_thresh)
        hop_thresh = float(hop_thresh)
        vel_adj = int(vel_adj)
        restart = int(restart)
        restart_freq = int(restart_freq)
        tera_mpi = int(tera_mpi)
        write_freq = int(write_freq)          
        print("Simulation parameters read from input.in:\n",
              "{} = {}\n".format("natoms", natoms),
              "{} = {}\n".format("maxsteps",maxsteps),
              "{} = {}\n".format("initial state", state),
              "{} = {}\n".format("timestep", dt),
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
    
    # ALLOCATE ALL NUMPY ARRAYS
    fx, fy, fz, fx_new, fy_new, fz_new, pot_eners, x_new, y_new, z_new = init_fep_arrays(natoms, nstates)
    
    # READ INITIAL/RESTART DATA
    if restart == 0:
        if init_vel:
            print("Initial velocities read from {}".format(vel_file_path)) 
        else: 
            print("No initial velocities.")
        print("Restart option turned OFF.")
             
        at_names, x, y, z  = read_geoms(natoms, geom_file_path)
        vx, vy, vz = read_velocs(init_vel, natoms, vel_file_path)

        # OBTAIN ATOMIC MASSES:
        am = assign_masses(at_names)
        
        # CENTER OF MASS REMOVAL 
        x, y, z = com_removal(x, y, z, am)
        
        # CALC INITIAL ENERGIES AND GRADIENTS
        if tera_mpi: 
            comm = tera_connect()  
            MO, CiVecs, NAC, blob, SMatrix, civec_size, nbf_size, blob_size, qmcharges, TDip, Dip = tera_init(comm, at_names, natoms, nstates, x, y, z)
             
        fx, fy, fz, pot_eners, \
        MO, CiVecs, blob = calc_forces(step, at_names, state, nstates, x, y, z, fx, fy, fz, pot_eners,
                                       ab_initio_file_path, tera_mpi, comm, sim_time, 
                                       MO, CiVecs, NAC, blob, SMatrix, civec_size, nbf_size, blob_size, qmcharges, TDip, Dip)

        pot_eners_array = np.copy(pot_eners)      
        
        Ekin, Epot, Etot, dE, dE_step = calc_energies(step, natoms, am,
                                                      state, pot_eners,
                                                      vx, vy, vz, Etot_init,
                                                      Etot_prev, ener_thresh)
        Etot_init = Etot
        init_step = 1
    else:
        rst_file_path = check_restart_files(restart, cwd)
        print(liner)
        init_step, at_names, state, \
        x, y, z, vx, vy, vz, fx, fy, fz, \
        Ekin, Epot, Etot, Etot_init, \
        pot_eners_array, CiVecs, MO, blob, civec_size, nbf_size, blob_size = read_restart(rst_file_path, natoms, nstates, tera_mpi)
         
        am = assign_masses(at_names)
        init_step = init_step + 1         #  main loop counter will start with a following step
        # if init_step = maxsteps:
        # error_exit(18, "Reaachem mxium number of steps")
        if tera_mpi:
            comm = tera_connect()  
            #  temp_ PREVENT overwrite of MO, CiVecs and blob data from restart 
            temp_MO, temp_CiVecs, NAC, temp_blob, SMatrix, civec_size, nbf_size, blob_size, qmcharges, TDip, Dip = tera_init(comm, at_names, natoms, nstates, x, y, z)
        
            # Check the shape of INIT and RESTART arrays - must equal
            # dim_match = "{}:{}\n{}:{}\n{}:{}\n ".format(temp_MO.shape, MO.shape, temp_CiVecs.shape, CiVecs.shape, temp_blob.shape, blob.shape)
            if not (temp_MO.shape == MO.shape and 
                    temp_CiVecs.shape == CiVecs.shape and
                    temp_blob.shape == blob.shape):
                error_exit(18, "dim_miss_match")
            # else:
            #    print("Restart and init dimensions match\n",dim_match)
                            
    check_output_file(cwd, natoms, restart, init_step, write_freq)

    print("{}\nInitial geometry:\nAt       X           Y            Z           MASS     [ANGSTROMS, AMU]".format(liner))
    xx = (x*BOHR_ANG).tolist()
    yy = (y*BOHR_ANG).tolist()
    zz = (z*BOHR_ANG).tolist()
    for iat in range(0, natoms):
        print(" {} {:12.8f}".format(at_names[iat], xx[iat]),
              "{:12.8f} {:12.8f}".format(yy[iat], zz[iat]),
              "{:12.8f}".format(am[iat]/AMU))

    print("Initial velocities:\nAt    VX      VY          VZ [ATOMIC UNITS]")
    for iat in range(0, natoms):
        print("".join("%2s" "   " "%6.4f" % (at_names[iat], vx[iat])),
              " %6.4f    %6.4f " % (vy[iat], vz[iat]))
   
    print("{}".format(liner),
          "\nSTEP    TIME/FS  DE_DRIFT/EV   DE_STEP/EV    HOP  STATE\n{}".format(liner)) 
    
    # open only once, otherwise intesive I/O disk trafic, especially for very fast ab initio calculations
    with open('movie.xyz', 'a') as mov_file, \
         open('energies.dat', 'a') as eners_file, \
         open('PES.dat', 'a') as pes_file, \
         open('velocities.dat', 'a') as vel_file, \
         open('state.dat', 'a') as state_file:

    #-------------------MAIN VERLET LOOP----------------------------------------------------------------
         for step in range(init_step, maxsteps + 1):
            sim_time = step * dt * AU_FS
            x_new, y_new, z_new = update_positions(dt, am, x, y, z, x_new, y_new, z_new, vx, vy, vz, fx, fy, fz)
    
            fx_new, fy_new, fz_new, pot_eners, \
            MO, CiVecs, blob = calc_forces(step, at_names, state, nstates, x_new, y_new, z_new, fx_new, fy_new, fz_new, pot_eners, ab_initio_file_path, 
                                           tera_mpi, comm, sim_time, MO, CiVecs, NAC, blob, SMatrix, civec_size, nbf_size, blob_size, qmcharges, TDip, Dip)
        
            if not method == "bomd":
                if step >= 2:
                    hop, outstate, v_scal_fac, prob = calc_hopp(method, state, pot_eners, pot_eners_array, Ekin, dt, hop_thresh)
    
                    if hop:
                        state = outstate
                        # use XYZ from prev. step to cacl F for a new state
                        
                        fx_new, fy_new, fz_new, pot_eners, \
                        MO, CiVecs, blob = calc_forces(step, at_names, state, nstates, x, y, z, fx_new, fy_new, fz_new, pot_eners, ab_initio_file_path,
                                                 tera_mpi,comm, sim_time, MO, CiVecs, NAC, blob, SMatrix, civec_size, nbf_size, blob_size, qmcharges, TDip, Dip)
    
                        #simple scaling or updatre velocities with new state forces
                        if not vel_adj:
                            vx, vy, vz = rescale_velocities(vx, vy, vz, v_scal_fac)
                        else:
                            vx, vy, vz = adjust_velocities(dt, am,
                                                           vx, vy, vz,
                                                           fx, fy, fz,
                                                           fx_new, fy_new, fz_new) 
    
                        # FXFYFZ for XYZ(t = hop) in the NEW state
                        fx = np.copy(fx_new)
                        fy = np.copy(fy_new)
                        fz = np.copy(fz_new)
                        EE = 0
                        for iat in range(0,natoms):
                            vvv = vx[iat] ** 2 + vy[iat] ** 2 + vz[iat] ** 2
                            EE = EE + (0.5 * am[iat] * vvv)
                        print("Old Ekin: {} \nScaled/Adjusted Ekin {}\n.".format(Ekin, EE))
    
                        # now finish the propagation step on new PES (FXFYFZ for new state)
                        x_new, y_new, z_new = update_positions(dt, am, x, y, z, x_new, y_new, z_new, 
                                                               vx, vy, vz, fx_new, fy_new, fz_new)
    
                        fx_new, fy_new, fz_new, pot_eners, \
                        MO, CiVecs, blob = calc_forces(step, at_names, state, nstates, x_new, y_new, z_new, fx_new, fy_new, fz_new, pot_eners, ab_initio_file_path,
                                                       tera_mpi,comm, sim_time, MO, CiVecs, NAC, blob, SMatrix, civec_size, nbf_size, blob_size, qmcharges, TDip, Dip)
                    #  Keep only last two pot_eners, 3rd pot_energs is appended during hopping evaluation for minima derivation, 
                    pot_eners_array = np.delete(pot_eners_array, 0, axis = 0) 
    
                #  Allways append the last pot_eners to the pot_eners_array 
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
    

            Etot_prev = Etot
            Ekin, Epot, Etot, dE, dE_step = calc_energies(step, natoms, am,
                                                          state, pot_eners,
                                                          vx, vy, vz, Etot_init,
                                                          Etot_prev, ener_thresh)
    
            # print("Ekin {}, Epot {}, Etot{}".format(Ekin, Epot, Etot))
            print(" {:<6d}  {:<7.4f}  {:<12.4f}".format(step, sim_time, dE * AU_EV),
                  " {:<13.4f}".format(dE_step * AU_EV),
                  "{}     {}".format(str(hop)[0], state))
            sys.stdout.flush()
            # SAVE POSITION, VELOCITIES, ENERGIES AND RESTART
            if (step%write_freq == 0):
                print_energies(step,write_freq, sim_time, Ekin, Epot, Etot, dE, dE_step, eners_file)
                print_pes(step, write_freq, sim_time, pot_eners, pes_file)
                print_positions(step, sim_time, natoms, at_names, x, y, z, mov_file)
                print_velocities(step, sim_time, natoms, at_names, vx, vy, vz, vel_file)
                print_state(step, write_freq, sim_time, state, state_file)
            
            if (step%restart_freq) == 0:
                print_restart(step, sim_time, natoms, at_names, state, timestep,
                              x, y, z, vx, vy, vz, fx, fy, fz, nstates,
                              Ekin, Epot, Etot, Etot_init, pot_eners_array, 
                              MO, CiVecs, blob, civec_size, nbf_size, blob_size, tera_mpi)

    # FINAL PRINTS
    print(liner)
    print("SIMULATION DONE")
    if tera_mpi:
        finish_tera(comm)   
    print("Output files: movie.xyz, velocities.xyz, PES.dat, energies.dat, state.dat")
    stopTime = datetime.now()
    calc_time = (datetime.now() - startTime)
    print("Simulation ended at: {}".format(stopTime))
    print("Overall simulation sim_time (hh:mm:ss): {}".format(calc_time))
    print(liner)
exit(0) 
