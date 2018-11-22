"""
Check input file: input.in, geom.in (input geometry), veloc.in pab initio program file,
Read input.in

each array (x,y,z) or (vx,vy,vz) or (fx,fy,fz)
holds NATOMS elements 
e.g. x[0] is x position of 1st atom, fx[1] force in x direction on 2 atom, etc.
"""

import sys
import os
import configparser
import numpy as np

try:
    from errors import error_exit
    from restart import (
        truncate_output_files, back_output_files
    )
    from constants import *
except ImportError as ime:
    print("Module {} not found.".format(ime.name))
    exit(1)

def file_check(cwd):
    # input files names - these are defaults otherwise not found
    veloc_file = "veloc.in"
    geom_file  = "geom.in"
    input_file = "input.in"
    
    # absolute path to the files  
    input_file_path = os.path.join(cwd, input_file)
    geom_file_path  = os.path.join(cwd, geom_file)
    vel_file_path = os.path.join(cwd, veloc_file)
    if (not os.path.isfile(input_file_path)):
         error_exit(0, " ")
    if (not os.path.isfile(geom_file_path)):
         error_exit(1, " ")
    if (not os.path.isfile(vel_file_path)):
        print("No initial velocities.")
        init_vel = 0
    else:
        print("Initial velocities read from {}".format(veloc_file)) 
        init_vel = 1
    
    return(input_file_path, geom_file_path, vel_file_path, init_vel)
    
def read_input(cwd, input_file_path):
    
    # Read parameter from input.in file
    cfg = configparser.ConfigParser(delimiters=('=', ':'),
        comment_prefixes=('#', ';'))
    cfg.read(input_file_path)       # Read file
    par=dict(cfg.items("Settings",raw=False))
    
    # To get rid of inline comments [0] and strip spaces
    print("Simulation will start with parameters: ")
    print("\n".join("{}: {}".format(k.strip(), v.split("#")[0].strip()) for k, v in par.items()))
    for p in par:
       par[p]=par[p].split("#",1)[0].strip(" ")  
        
    # Now that we know everything - check for ab initio interface
    if not par['tera_mpi']:
        abinit_file = "ABINITIO/{}".format(par['abinitio'])
        ab_initio_file_path  = os.path.join(cwd, abinit_file)
    else 
        abinit_file = "ABINITIO/{}".format(par['abinitio'])
        ab_initio_file_path  = os.path.join(cwd, abinit_file)
    if (not os.path.isfile(ab_initio_file_path)):
         error_exit(5, " ")
         
    return(par, ab_initio_file_path)

def check_output_file(cwd, natoms, restart, init_step):
    if (restart == 1):
        back_output_files()
    elif (restart > 1):
        backup_output_files()
        truncate_output_files(init_step, natoms):
    elif (restart == 0):
        if (os.path.isfile("movie.xyz")):
         error_exit(8, "movie.xyz")
        if (os.path.isfile("velocities.xyz")):
         error_exit(8, "velocities.xyz")     
        if (os.path.isfile("energies.dat")):
         error_exit(8, "energies.dat")  
        if (os.path.isfile("PES.dat")):
         error_exit(8, "PES.dat")     
        if (os.path.isfile("state.dat")):
         error_exit(8, "state.dat")
        if (os.path.isfile("restart.in")):
         print("File restart.in exists, but restart option is turned off.")
         error_exit(8, "restart.in")
    try:
        with open('movie.xyz', 'a') as mov_file, \
             open('energies.dat', 'a') as eners_file, \
             open('PES.dat', 'a') as pes_file, \
             open('velocities.dat', 'a') as vel_file, \
             open('state.dat', 'a') as state_file:
    except IOError as ioe:
         print("{Failed to open output file})".format(ioe.strerror))
         # TODO MPI exit
    return(mov_file, eners_file, pes_file, vel_file, state_file)   
def output_files_close(mov_file, eners_file, pes_file, vel_file, state_file):
    mov_file.closed
    eners_file.closed
    pes_file.closed
    vel_file.closed
    state_file.closed
    return
def read_geoms(natoms, geom_file_path):
    #if restart == 1 : read  last two geoms
    x = np.zeros(natoms, dtype=np.float64)  
    y = np.zeros(natoms, dtype=np.float64)  
    z = np.zeros(natoms, dtype=np.float64)  

    at_names = []
 # READ INITIAL POSITIONS:   
    with open(geom_file_path, 'r') as igf:  # igf input geom file 
     atoms = igf.readline()  # first line in geom file is number of atoms
     if not (int(atoms) == natoms):  error_exit(2, " ")                         
     garbage = igf.readline()  # comment 
     for iat in range(0, natoms):
        line = igf.readline().split()
        at_names.append(capitalize_2th(str(line[0])))
        x[iat] = np.float64(line[1]) * ang_bohr  # atomic units : Bohr
        y[iat] = np.float64(line[2]) * ang_bohr
        z[iat] = np.float64(line[3]) * ang_bohr
     
     igf.close()
     return(at_names, x, y, z)

 # READ INITIAL VELOCITIES:
def read_velocs(init_vel, natoms, vel_file_path):  
    vx = np.zeros(natoms, dtype=np.float64)  
    vy = np.zeros(natoms, dtype=np.float64)  
    vz = np.zeros(natoms, dtype=np.float64)  
    if init_vel:
        with open(vel_file_path, 'r') as ivf:  
            atoms = ivf.readline()  # first line is number of atoms
            if not (int(atoms) == natoms):  
                error_exit(2, " ")                         
            garbage = ivf.readline()  # second comment line
            for iat in range(0, natoms):
                line = ivf.readline().split()
                vx[iat] = float(line[1])
                vy[iat] = float(line[2])
                vz[iat] = float(line[3])
        ivf.close()
    return(vx, vy, vz)

def capitalize_2th(s):
    # Capitalizace first letter, lower second -
    # avoid problems with different name for atoms in ab initio codes
    return s[:1].capitalize() + s[1:].lower()

def init_fep_arrays(natoms, nstates):
    # Initialize empty forces, energies and pos_new arrays
    
    #f(t) f_new(t+dt)
    fx = np.zeros(natoms, dtype = np.float64)  
    fy = np.zeros(natoms, dtype = np.float64)   
    fz = np.zeros(natoms, dtype = np.float64)  
    fx_new = np.zeros(natoms, dtype = np.float64)    
    fy_new = np.zeros(natoms, dtype = np.float64)  
    fz_new = np.zeros(natoms, dtype = np.float64)  

    # potential energie array
    pot_eners = np.zeros(nstates, dtype = np.float64)   
    
    # new positions in shift_pos in verlet step
    x_new = np.zeros(natoms, dtype=np.float64)  
    y_new = np.zeros(natoms, dtype=np.float64)  
    z_new = np.zeros(natoms, dtype=np.float64)  
        
    return(fx, fy, fz, fx_new, fy_new, fz_new, pot_eners, x_new, y_new, z_new)
    

    
def com_removal(x, y, z, am):
    totmass, xsum, ysum, zsum = 0.0, 0.0, 0.0, 0.0
    for iat in range(0,len(x)):
        xsum += x[iat] * am[iat] 
        ysum += y[iat] * am[iat] 
        zsum += z[iat] * am[iat] 
        totmass += am[iat] 
        
    xcom = xsum / totmass
    ycom = ysum / totmass
    zcom = zsum / totmass
     
    for iat in range(0,len(x)):
        x[iat] = x[iat] - xcom
        y[iat] = y[iat] - ycom
        z[iat] = z[iat] - zcom 
    
    return(x,y,z)
