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
from errors import error_exit


au_fs = 0.02418884326505e0      #atomic units to femtosecs
au_eV = 27.21139
amu   = 1822.888484264545e0             # atomic mass unit  me = 1 AMU*atomic weight
ang_bohr = 1.889726132873e0     # agstroms to bohrs

def file_check(cwd):
    # input files names - these are defaults otherwise not found
    veloc_file = "veloc.in"
    geom_file  = "geom.in"
    input_file = "input.in"
    
    # absolute path to the files  
    input_file_path = os.path.join(cwd,input_file)
    geom_file_path  = os.path.join(cwd,geom_file)
    vel_file_path = os.path.join(cwd,veloc_file)
    if (not os.path.isfile(input_file_path)):
         error_exit(0)
    if (not os.path.isfile(geom_file_path)):
         error_exit(1)
    if (not os.path.isfile(vel_file_path)):
        print("No initial velocities.")
        init_vel = 0
    else:
        print("Initial velocities read from {}".format(veloc_file)) 
        init_vel = 1
    
    return(input_file_path, geom_file_path, vel_file_path, init_vel)
    
def read_input(cwd,input_file_path):
    
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
    abinit_file = "ABINITIO/{}.sh".format(par['abinitio'])
    ab_initio_file_path  = os.path.join(cwd,abinit_file)
    
    if (not os.path.isfile(ab_initio_file_path)):
         error_exit(5)
         
    return(par,ab_initio_file_path)

def read_geoms(natoms,geom_file_path):
    #if restart == 1 : read  last two geoms
    # could be xyz matrix [natoms,3] as well, but this should be more clear reading
    x = np.zeros(natoms,dtype=np.float64)  
    y = np.zeros(natoms,dtype=np.float64)  
    z = np.zeros(natoms,dtype=np.float64)  
    at_names = []
 # READ INITIAL POSITIONS:   
    with open(geom_file_path,'r') as igf:  # igf input geom file 
     atoms = igf.readline()  # first line in geom file is number of atoms
     if not (int(atoms) == natoms):  error_exit(2)                         
     garbage = igf.readline()  # comment 
     for iat in range(0,natoms):
        line = igf.readline().split()
        at_names.append(capitalize_2th(str(line[0])))
        x[iat] = np.float64(line[1]) * ang_bohr  # atomic units : Bohr
        y[iat] = np.float64(line[2]) * ang_bohr
        z[iat] = np.float64(line[3]) * ang_bohr
     
     igf.close()
     return(at_names,x,y,z)
 # READ INITIAL VELOCITIES:
def read_velocs(init_vel, natoms, vel_file_path):  
    vx = np.zeros(natoms,dtype=np.float64)  
    vy = np.zeros(natoms,dtype=np.float64)  
    vz = np.zeros(natoms,dtype=np.float64)  
    if init_vel:
        with open(vel_file_path,'r') as ivf:  
            atoms = ivf.readline()  # first line is number of atoms
            if not (int(atoms) == natoms):  
                error_exit(2)                         
            garbage = ivf.readline()  # second comment line
            for iat in range(0,natoms):
                line = ivf.readline().split()
                vx[iat] = float(line[1])
                vy[iat] = float(line[2])
                vz[iat] = float(line[3])
        ivf.close()
    return(vx,vy,vz)

# Capitalizace first letter, lower second - avoid problems with different name for atoms
def capitalize_2th(s):
    return s[:1].capitalize() + s[1:].lower()

def init_forces_energs(natoms,nstates):
# Initialize empty forces array
    
    #f(t)
    fx = np.zeros(natoms,dtype=np.float64)  
    fy = np.zeros(natoms,dtype=np.float64)   
    fz = np.zeros(natoms,dtype=np.float64)  
    #f_new(t+dt)
    fx_new = np.zeros(natoms,dtype=np.float64)    
    fy_new = np.zeros(natoms,dtype=np.float64)  
    fz_new = np.zeros(natoms,dtype=np.float64)  
    pot_eners = [0.00000000 ] * nstates # potential energy from ab initio calculations
    
    return(fx, fy, fz, fx_new, fy_new, fz_new, pot_eners)
    
def com_removal(x,y,z,am):
    totmass, xsum, ysum, zsum = 0.0, 0.0, 0.0, 0.0
    for iat in range(0,len(x)):
        xsum += x[iat] * am[iat] 
        ysum += y[iat] * am[iat] 
        zsum += z[iat] * am[iat] 
        totmass += am[iat] 
        
    xcom = xsum / totmass
    ycom = xsum / totmass
    zcom = xsum / totmass
     
    for iat in range(0,len(x)):
        x[iat] = x[iat] - xcom
        y[iat] = y[iat] - ycom
        z[iat] = z[iat] - zcom 
    
    return(x,y,z)

