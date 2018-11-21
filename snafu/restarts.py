import numpy as np
import os
import shutil
import re
import subprocess

try:
    from errors import error_exit
    from constants import *
except ImportError as ime:
    print("Module {} not found.".format(ime.name))
    exit(1)

def check_restart_files(restart, cwd):
    rst_file = "restart.in"
    rst_file_path = os.path.join(cwd, rst_file)
    if restart < 0:
            error_exit(11, "(restart < 0)")
    #  elif restart == 0 dealt in check_output_files
    elif restart == 1:
        rst_file = "restart.in"  
        rst_file_path = os.path.join(cwd, rst_file)
        if (os.path.isfile(rst_file_path)):
            print("Restart={}, start from the last completed".format(restart),
                  "step.\nRestart file {} found.".format(rst_file))
        else:
            error_exit(10, "({})".format(rst_file))
    elif restart > 1:
        rst_file = "restart_{}.in".format(restart)
        rst_file_path = os.path.join(cwd, rst_file)
        if (os.path.isfile(rst_file_path)):
            print("Restart={}, start from {} step".format(restart, restart),
                  "Restart file {} found".format(rst_file))
        else:
            error_exit(10, "({})".format(rst_file))
    return(rst_file_path)
    
def read_restart(rst_file_path, natoms):
    #  with open restart.in as rsf:
    #  print(rst_file_path)
    rst_file = rst_file_path
    with open(rst_file_path, 'r') as rstf:
        for num, line in enumerate(rstf):
            print(num, line)        
            if re.search(r'Step', line):
                step = int(line.split()[1])
            if re.search(r'State', line):
                state = int(line.split()[1])
            if re.search(r'Ekin', line):
                Ekin = float(line.split()[1])
            if re.search(r'Epot', line):
                Epot = float(line.split()[1])
            if re.search(r'Etot', line):
                Etot = float(line.split()[1])
            if re.search(r'Etot_init', line):
               Etot_init = float(line.split()[1])
            if re.search(r'Positions', line):
               pnum = num
            if re.search(r'Velocities', line):
                vnum = num
            if re.search(r'Forces', line):
                fnum = num
            if re.search(r'Pot_eners_array', line):
                peanum = num    
        
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
                   
    rstf.closed
    #pot_eners_array = np.loadtxt(rst_file, dtype=np.float64, delimiter=None, skiprows=8)
    np.set_printoptions(precision=10, formatter={'float': '{: 0.8f}'.format})        
    print(x)
    print(y)
    print(z)
    print(vx)
    print(vy)
    print(vz)
    print(fx)
    print(fy)
    print(fz)
    print(at_names)
    print(pot_eners_array)
    return(step, at_names, state, 
           x, y, z, vx, vy, vz, 
           fx, fy, fz, Ekin, Epot, Etot, Etot_init, pot_eners_array)

def print_restart(
        step, time, natoms, at_names, state, timestep,
        x, y, z, vx, vy, vz, fx, fy, fz,
        Ekin, Epot, Etot, Etot_init, pot_eners_array,
        restart_freq):

    inits_line = ("Step: {:d}".format(step),
                "State: {:d}".format(state),
                "Natoms: {:d}".format(natoms),
                "Ekin: {:14.10f}".format(Ekin),
                "Epot: {:14.10f}".format(Epot),
                "Etot: {:14.10f}".format(Etot),
                "Etot_init: {:14.10f}".format(Etot_init),
                "Pot_eners_array:\n"
                )

    rst_file = "restart.in"
    rsf.truncate(0)  #  empty file
    
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
    rsf.closed

    if not (step%restart_freq):
        rst_file = "restart_{}.in".format(step)
        print("Writing restart information to {} file.".format(rst_file))
        shutil.copy("restart.in", rst_file, follow_symlinks=True) 
      
    return()
    
def truncate_output_files(init_step, natoms):
# movie.xyz energies.dat input.in state.dat snafu.out velocities.xyz PES.dat 
# init_step is read from restart file
    print("Restart from {} step. Output file will be truncated".format(init_step))
    natoms_lines = (natoms+2)*init_step   #  header
    step_lines = init_step                
    input_files = []
    input_files.append(["movie.xyz", natoms_lines])
    input_files.append(["velocities.xyz", natoms_lines])
    input_files.append(["PES.dat", step_lines])
    input_files.append(["energies.dat", step_lines])
    input_files.append(["state.dat", step_lines])
    
    for xx in range(len(input_files)):
        nlines = input_files[xx][1]
        input_file = input_files[xx][0]
        cmd_trunc = "head -n{} {} > temp_file; mv temp_file {} ".format(nlines, input_file, input_file)    
        try:
            trim_file = subprocess.run(cmd_trunc, stdout=None, stderr=subprocess.PIPE, shell = True, check = True)	
        except subprocess.CalledProcessError as cpe: 
            print("Return code: {}\nError: {}".format(cpe.returncode, cpe.stderr))
            exit(1)
        else:
            print("File {} was truncated after {} steps.".format(input_file,init_step))
    return()

def backup_output_files(cwd):
    # RESTART PART - CREATING BACKUP OF OUTPUT FILES
        
    N=0
    while os.path.isdir(os.path.join(cwd,"PREV_RUN"+str(N))):
        print("{} backup folder already exists".format("PREV_RUN"+str(N)))
        N += 1
    else:
        backup_folder = os.path.join(cwd,"PREV_RUN"+str(N))
        os.mkdir(backup_folder)
        print("Creating backup folder {}".format("PREV_RUN"+str(N)))
    output_files= []
    output_files.append("movie.xyz")
    output_files.append("velocities.xyz")
    output_files.append("PES.dat")
    output_files.append("energies.dat")
    output_files.append("state.dat")
    output_files.append("snafu.out")
    output_files.append("restart.in")
    restart_files = [ rf for rf in rest_files if re.search(r'restart', rf)]
    backup_files = output_files + restart_files
    for bf in backup_files:
        shutil.copy(bf,backup_folder)
    print("Output data were backed-up.\n"format(os.listdir(backup_folder)))
    return()