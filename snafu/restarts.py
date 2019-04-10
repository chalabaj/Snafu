import numpy as np
import os
import shutil
import re
import subprocess
import sys
current_module = sys.modules[__name__]
from defaults import liner
try:
    from errors import error_exit
    from constants import *  
except ImportError as ime:
    error_exit(19, "Module {} in {} not found.".format(ime,current_module))
try:
    tera_mpi = int(os.environ['MPI_TERA'])
except KeyError as ke:
    print("MPI_TERA variable was not exported, assuming MPI_TERA=0. Warning: this may cause deadlock if MPI has been already initiated")
    tera_mpi = 0
if tera_mpi:
    from tera_propagates import (global_except_hook)
    sys.excepthook = global_except_hook  


# ---------------------------------------------------------------------------------------------------------------------------------------------         
def check_restart_files(restart, cwd):
    rst_file = "restart.in"
    rst_file_path = os.path.join(cwd, rst_file)
    if restart < 0:
        error_exit(11, "(restart < 0)")
    #  elif restart == 0 dealt in check_output_files
    #elif restart == 1:
    #    rst_file_path = os.path.join(cwd, rst_file)
    #    if (os.path.isfile(rst_file_path)):
    #        print("Restart={}, start from the last completed".format(restart),
    #              "step.\nRestart file {} found.".format(rst_file))
    #    else:
    #        error_exit(10, "({})".format(rst_file))
    elif restart > 0:
        rst_file = "restart_{}.in".format(restart)
        rst_file_path = os.path.join(cwd, rst_file)
        if (os.path.isfile(rst_file_path)):
            print("Restart={}".format(restart),
                  "\nRestart file {} found".format(rst_file))
        else:
            error_exit(10, "({})".format(rst_file))
    return(rst_file_path)
    
def read_restart(rst_file_path, natoms, nstates, tera_mpi):
    # TO DO: should be done with far less condiotion, however there would have to be another convertion from disctionary to variables of different types
    err = "0"
    rst_file = rst_file_path
    #step, state, Ekin, Epot, Etot, Etot_init, pnum, vnum, fnum, peanum, blob_size, civec_size, nbf_size, monum, civecnum, blobnum = None
    with open(rst_file_path, 'r') as rstf:
        for num, line in enumerate(rstf):
            if re.search(r'Step', str(line)):
                step = int(line.split()[1])
            if re.search(r'State', str(line)):
                state = int(line.split()[1])
            if re.search(r'Ekin', str(line)):
                Ekin = float(line.split()[1])
            if re.search(r'Epot', str(line)):
                Epot = float(line.split()[1])
            if re.search(r'Etot', str(line)):
                Etot = float(line.split()[1])
            if re.search(r'Etot_init', str(line)):
                Etot_init = float(line.split()[1])   
            if re.search(r'Positions', str(line)):
                pnum = num
            if re.search(r'Velocities', str(line)):
                vnum = num
            if re.search(r'Forces', str(line)):
                fnum = num
            if re.search(r'Pot_eners_array', str(line)):
                peanum = num 
            # TC WF
            if re.search(r'blob_size', str(line)):
                blob_size= int(line.split()[1])
            if re.search(r'civec_size', str(line)):
                civec_size = int(line.split()[1])  
            if re.search(r'nbf_size', str(line)):
                nbf_size= int(line.split()[1])  
            if re.search(r'MO', str(line)):
                monum = num
            if re.search(r'CiVecs', str(line)):
                civecnum = num
            if (re.search(r'Blob:$', str(line))):
                blobnum = num
        try: 
            if not None in (step, state, Ekin, Epot, Etot, Etot_init, pnum, vnum, fnum, peanum, blob_size, civec_size, nbf_size, monum, civecnum, blobnum):
                print("ALL RESTARTION OPTIONS FOUND. READING RESTART DATA:")  
        except Exception as ne:
            error_exit(17, str(ne))
        try:          
            atnames = np.genfromtxt(rst_file, dtype=np.dtype('str'),
                                    skip_header=pnum+1, max_rows=natoms, usecols=[0])
            at_names =  atnames.tolist()
            
            pot_eners_array = np.genfromtxt(rst_file, dtype=np.float64,
                                            skip_header=peanum+1, max_rows=2,
                                            usecols=[0,1,2])         
            fx = np.genfromtxt(rst_file, dtype=np.float64,
                               skip_header=fnum+1, max_rows=natoms, usecols=[1])  
            fy = np.genfromtxt(rst_file, dtype=np.float64,
                               skip_header=fnum+1, max_rows=natoms, usecols=[2])
            fz = np.genfromtxt(rst_file, dtype=np.float64,
                               skip_header=fnum+1, max_rows=natoms, usecols=[3])
                                       
                                       
            x = np.genfromtxt(rst_file, dtype=np.float64,
                                      skip_header=pnum+1, max_rows=natoms, usecols=[1])  
            y = np.genfromtxt(rst_file, dtype=np.float64, 
                                      skip_header=pnum+1, max_rows=natoms, usecols=[2])
            z = np.genfromtxt(rst_file, dtype=np.float64,
                                      skip_header=pnum+1, max_rows=natoms, usecols=[3])
                                      
            vx = np.genfromtxt(rst_file, dtype=np.float64,
                                       skip_header=vnum+1, max_rows=natoms, usecols=[1])  
            vy = np.genfromtxt(rst_file, dtype=np.float64,
                                       skip_header=vnum+1, max_rows=natoms, usecols=[2])
            vz = np.genfromtxt(rst_file, dtype=np.float64,
                                       skip_header=vnum+1, max_rows=natoms, usecols=[3])
                                       
            if not tera_mpi == 0:
                #CiVecs = np.zeros((civec_size, nstates),dtype=np.float64)                   
                CiVecs = np.genfromtxt(rst_file, dtype=np.float64,
                                           skip_header=civecnum+1, max_rows=civec_size, usecols=(x for x in range(nstates)))
                #MO = np.zeros((nbf_size, nbf_size),dtype=np.float64)
                MO = np.genfromtxt(rst_file, dtype=np.float64,
                                           skip_header=monum+1, max_rows=nbf_size, usecols=(x for x in range(nbf_size))) 
                #blob = np.zeros((blob_size),dtype=np.float64)
                blob = np.genfromtxt(rst_file, dtype=np.float64,
                                           skip_header=blobnum+1, max_rows=blob_size, usecols=0)                             
        except Exception as expt:
            print(expt)
            error_exit(16)
           
    rstf.closed
    #pot_eners_array = np.loadtxt(rst_file, dtype=np.float64, delimiter=None, skiprows=8)
    np.set_printoptions(precision=10, formatter={'float': '{: 0.8f}'.format})        
    print("Reading restart data:")
    print(at_names)
    print("XYZ:\n{}".format(np.dstack((x,y,z))))
    print("VX VY VZ:\n{}".format(np.dstack((vx,vy,vz))))
    print("FX FY FZ:\n{}".format(np.dstack((fx,fy,fz))))
    print("Pot ener step-1/step \n{}".format(pot_eners_array))
    return(step, at_names, state, x, y, z, vx, vy, vz, 
           fx, fy, fz, Ekin, Epot, Etot, Etot_init, pot_eners_array, CiVecs, MO, blob, civec_size, nbf_size, blob_size)

def print_restart(step, sim_time, natoms, at_names, state, timestep,
                  x, y, z, vx, vy, vz, fx, fy, fz, nstates,
                  Ekin, Epot, Etot, Etot_init, pot_eners_array, 
                  MO, CiVecs, blob, civec_size, nbf_size, blob_size, tera_mpi):
    
    rst_file_step = "restart_{}.in".format(step)
    #print("Writing restart information to {} file.".format(rst_file_step))
    with open(rst_file_step, 'w') as rsf:
        # if it existed, empty it
        # rsf.truncate(0)  #  empty file
        # rsf.seek(0)
        inits_line = ("Step: {:d}".format(step),
                    "State: {:d}".format(state),
                    "Natoms: {:d}".format(natoms),
                    "Ekin: {:14.10f}".format(Ekin),
                    "Epot: {:14.10f}".format(Epot),
                    "Etot: {:14.10f}".format(Etot),
                    "Etot_init: {:14.10f}".format(Etot_init),
                    "civec_size: {}".format(int(civec_size)), 
                    "blob_size: {}".format(int(blob_size)),
                    "nbf_size: {}".format(int(nbf_size)), 
                    "Pot_eners_array:\n"
                    )
        
        rsf.write('\n'.join(inits_line))
    
        np.savetxt(rsf, pot_eners_array, fmt="%20.10f", delimiter=' ', newline='\n')
    
        rsf.write("Positions AT X Y Z (Bohrs):\n")
        xx = x.tolist()
        yy = y.tolist()
        zz = z.tolist()
        for iat in range(0, natoms):     
            p_line = "{} {:20.10f} {:20.10f} {:20.10f}\n".format(at_names[iat],
                                                                 xx[iat],
                                                                 yy[iat],
                                                                 zz[iat])
            rsf.write(p_line)
                
        rsf.write("Velocities: AT VX VY VZ (a.u.):\n")
        vvx = vx.tolist()
        vvy = vy.tolist()
        vvz = vz.tolist()
        for iat in range(0, natoms):     
            v_line = "{} {:20.10f} {:20.10f} {:20.10f}\n".format(at_names[iat],
                                                                 vvx[iat],
                                                                 vvy[iat],
                                                                 vvz[iat])
            rsf.write(v_line)
    
        rsf.write("Forces: AT FX FY FZ (a.u.):\n")
        ffx = fx.tolist()
        ffy = fy.tolist()
        ffz = fz.tolist()
        for iat in range(0, natoms):
                f_line = "{} {:20.10f} {:20.10f} {:20.10f}\n".format(at_names[iat],
                                                                   ffx[iat],
                                                                   ffy[iat],
                                                                   ffz[iat])
                rsf.write(f_line)
                
        if not tera_mpi == 0:
            rsf.write("MO:\n")
            np.savetxt(rsf, MO, fmt="%20.10f", delimiter=' ', newline='\n')
            
            rsf.write("CiVecs:\n")
            np.savetxt(rsf, CiVecs, fmt="%20.10f", delimiter=' ', newline='\n')
            
            rsf.write("Blob:\n")
            np.savetxt(rsf, blob, fmt="%20.10f", delimiter=' ', newline='\n')
    
        rsf.flush()  # OLD: need to flush, cause it is still open and buffer probably not full yet
        
    return()
    
def truncate_output_files(init_step, write_freq, natoms):
# movie.xyz energies.dat input.in state.dat snafu.out velocities.xyz PES.dat 
# init_step is read from restart file
    print("\nOutput files will be truncated".format(init_step))
    natoms_lines = (natoms+2)*init_step/write_freq  #  header
    step_lines = (init_step/write_freq)+1     
           
    input_files = []
    input_files.append(["movie.xyz", natoms_lines])
    input_files.append(["velocities.dat", natoms_lines])
    input_files.append(["PES.dat", step_lines])
    input_files.append(["energies.dat", step_lines])
    input_files.append(["state.dat", step_lines])
    
    for xx in range(len(input_files)):
        nlines = int(input_files[xx][1])
        input_file = input_files[xx][0]
        cmd_trunc = "head -n{} {} > temp_file && mv temp_file {} ".format(nlines, input_file, input_file)  #  && = and if True
        try:
            trim_file = subprocess.run(cmd_trunc, stdout=None, stderr=subprocess.PIPE, shell = True, check = True)	
        except subprocess.CalledProcessError as cpe: 
            print("Warning: error during output {} file truncations process: \n{}".format(input_file,cpe.stderr),
                  "\nFile was probably moved somewhere else.")
        else:
            print("File {} was truncated after {} steps.".format(input_file,init_step))
    return()

def backup_output_files(cwd):
    # RESTART PART - CREATING BACKUP OF OUTPUT FILES
    print(liner)    
    N=0
    while os.path.isdir(os.path.join(cwd,"PREV_RUN"+str(N))):
        print("{} backup folder already exists".format("PREV_RUN"+str(N)))
        N += 1
    else:
        backup_folder = os.path.join(cwd,"PREV_RUN"+str(N))
        os.mkdir(backup_folder)
        print("Creating backup folder: {}\n".format("PREV_RUN"+str(N)))
    output_files= []
    output_files.append("movie.xyz")
    output_files.append("velocities.dat")
    output_files.append("PES.dat")
    output_files.append("energies.dat")
    output_files.append("state.dat")
    #  output_files.append("snafu.out") file is rewritten upon restart
    output_files.append("restart.in")
    rest_files = os.listdir(os.getcwd())
    restart_files = [ rf for rf in rest_files if re.search(r'restart', rf)]
    backup_files = output_files + restart_files
    for bf in backup_files:
        try:
            shutil.copy(bf,backup_folder)
            print("File {} was not backed-up.".format(bf))
        except FileNotFoundError as FNT:
            print("Warning: file {} was not found and was not be backed-up!!!".format(FNT.filename))
    print("Old output data were backed-up:\n{}".format(os.listdir(backup_folder)))
    return()
