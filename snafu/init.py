"""
Check input file: input.in, geom.in (input geometry), veloc.in pab initio program file,
Read input.in
"""
import numpy as np
import sys
import os
import configparser
from snafu.errors import error_exit

def file_check(cwd):
    
    veloc_file = "veloc.in"
    geom_file  = "geom.in"
    input_file = "input.in"
       
    input_file_path = os.path.join(cwd,input_file)
    geom_file_path  = os.path.join(cwd,geom_file)
    veloc_file_path = os.path.join(cwd,veloc_file)
    
    if (not os.path.isfile(input_file_path)):
         error_exit(0)
    if (not os.path.isfile(geom_file_path)):
         error_exit(1)
    if(not os.path.isfile(veloc_file_path)):
         print("No initial velocities, continue with zero initial velocities.")
         veloc_init = 0
    else: veloc_init = 1
    
    return(input_file_path, geom_file_path, veloc_file_path, veloc_init)
    
def read_input(input_file_path):
    cfg = configparser.ConfigParser(delimiters=('=', ':'),comment_prefixes=('#', ';'))
    cfg.read(input_file_path)       # Read file
    par=dict(cfg.items("Settings",raw=False))
    
    # To get rid of inline comments [0] and strip spaces
    print("Simulation will start with parameters: ")
    print("\n".join("{}: {}".format(k.strip(), v.split("#")[0].strip()) for k, v in par.items()))
    for p in par:
       par[p]=par[p].split("#",1)[0].strip(" ")  
    return(par)

def read_geoms(natoms,geom_file_path):
    #if restart == 1 : read  last two geoms
    # could be xyz matrix [natoms,3] as well, but this should be more clear reading
    x = np.zeros(shape=(natoms,1)) 
    y = np.zeros(shape=(natoms,1)) 
    z = np.zeros(shape=(natoms,1)) 
    at_names = []
 # READ INITIAL POSITIONS:   
    with open(geom_file_path,'r') as igf:  # igf input geom file 
     atoms = igf.readline()  # first line in geom file is number of atoms
     if not (int(atoms) == natoms):  error_exit(2)                         
     garbage = igf.readline()  # comment 
     for iat in range(0,natoms):
        line = igf.readline().split()
        at_names.append(capitalize_2th(str(line[0])))
        x[iat] = float(line[1])
        y[iat] = float(line[2])
        z[iat] = float(line[3])
     
     igf.close()
     return(at_names,x,y,z)
 # READ INITIAL VELOCITIES:
def read_velocs(veloc_init,natoms,veloc_file_path):  
    vx = np.zeros(shape=(natoms,1)) 
    vy = np.zeros(shape=(natoms,1)) 
    vz = np.zeros(shape=(natoms,1)) 
    if veloc_init == 1:
     with open(veloc_file_path,'r') as ivf:  # ivf input veloc file 
       atoms = ivf.readline()  # first line in geom file is number of atoms
       if not (int(atoms) == natoms):  error_exit(2)                         
       garbage = ivf.readline()  # comment in file
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