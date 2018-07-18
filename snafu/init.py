"""
Check input file: input.in, geom.in (input geometry), veloc.in pab initio program file,
Read input.in
"""

import sys
import os
import configparser
from snafu.errors import error_exit

au_fs = 0.02418884326505      #atomic units to femtosecs
au_eV = 27.21139
amu   = 1822.8885             # atomic mass unit  me = 1 AMU*atomic weight
ang_bohr = 1.889726132873     # agstroms to bohrs

def file_check(cwd):
    # input files names - these are defaults otherwise not found
    veloc_file = "veloc.in"
    geom_file  = "geom.in"
    input_file = "input.in"
    
    # absolute path to the files  
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
    
def read_input(cwd,input_file_path):
    
    # Read parameter from input.in file
    cfg = configparser.ConfigParser(delimiters=('=', ':'),comment_prefixes=('#', ';'))
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
    x = [ 0.0 ] * natoms 
    y = [ 0.0 ] * natoms 
    z = [ 0.0 ] * natoms
    at_names = []
 # READ INITIAL POSITIONS:   
    with open(geom_file_path,'r') as igf:  # igf input geom file 
     atoms = igf.readline()  # first line in geom file is number of atoms
     if not (int(atoms) == natoms):  error_exit(2)                         
     garbage = igf.readline()  # comment 
     for iat in range(0,natoms):
        line = igf.readline().split()
        at_names.append(capitalize_2th(str(line[0])))
        x[iat] = float(line[1]) * ang_bohr  # atomic units : Bohr
        y[iat] = float(line[2]) * ang_bohr
        z[iat] = float(line[3]) * ang_bohr
     
     igf.close()
     return(at_names,x,y,z)
 # READ INITIAL VELOCITIES:
def read_velocs(veloc_init,natoms,veloc_file_path):  
    vx = [ 0.0 ] * natoms
    vy = [ 0.0 ] * natoms 
    vz = [ 0.0 ] * natoms 
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

def create_output_file(files):

    try:
        for x in files:
          cmd = "touch {}".format(x)
          os.system(cmd)
    except OSError:
       error_exit(3)
    except WindowsError:
       error_exit(3)
    return()   

def init_forces_potenergs(natoms,nstates):
# Initialize empty forces array
    
    #f(t)
    fx =([ 0.0 ] * natoms)  
    fy =([ 0.0 ] * natoms)  
    fz =([ 0.0 ] * natoms)
    #f_new(t+dt)
    fx_new = [ 0.0 ] * natoms  
    fy_new = [ 0.0 ] * natoms
    fz_new = [ 0.0 ] * natoms
    pot_eners = [0.0] * nstates # potential energy from ab initio calculations
    
    return(fx,fy,fz,fx_new,fy_new,fz_new, pot_eners)
