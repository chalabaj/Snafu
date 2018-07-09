"""
Check input file: input.in, geom.in (input geometry), veloc.in pab initio program file,
Read input.in
"""

import sys
import os
import configparser

def file_check(cwd):
    print("Working directory: ",cwd)
    veloc_file = "veloc.in"
    geom_file  = "geom.in"
    input_file = "input.in"
    
    print("Cheking input files.\n")
    
    input_file_path = os.path.join(cwd,input_file)
    geom_file_path  = os.path.join(cwd,geom_file)
    veloc_file_path = os.path.join(cwd,veloc_file)
    
    if (not os.path.isfile(input_file_path)) or (not os.path.isfile(geom_file_path)):
         print("File input.in or geom.in does not exist. Exiting... \n")
         sys.exit(1)
    if(not os.path.isfile(veloc_file_path)):
         print("File with initial velocities not found, continue with zero initial velocities.")
         veloc_init = 0
    else: veloc_init = 1
    print("Files check ok.\n")
    return(input_file_path, geom_file_path, veloc_file_path, veloc_init)
    
def read_input(input_file_path):
    cfg = configparser.ConfigParser(delimiters=('=', ':'),comment_prefixes=('#', ';'))
    cfg.read(input_file_path)       # Read file
    par=dict(cfg.items("Settings",raw=False))
    
    # To get rid of inline comments [0] and strip spaces
    print("Starting simulation with folowing parameters: ")
    print("\n".join("{}: {}".format(k.strip(), v.split("#")[0].strip()) for k, v in par.items()))
    for p in par:
       par[p]=par[p].split("#",1)[0].strip(" ")  
    return(par)