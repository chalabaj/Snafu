"""
Velocity verlet integrator scheme
"""

import math
import sys
import os
import random
import time
import subprocess
import re
import numpy as np

try:
    from errors import error_exit
    from prints import print_energies, print_pes
    from constants import *
except ImportError as ime:
    print("Module {} not found.".format(ime.name))
    exit(1)
        
# CONSTANTS
au_fs = 0.02418884326505e0   # atomic units to femtosecs
au_eV = 27.21139
amu = 1822.888484264545e0    # atomic mass unit me = 1 AMU*atomic weight
ang_bohr = 1.889726132873e0  # agstroms to bohrs

def update_positions(dt, am, x, y, z, vx, vy, vz, fx, fy, fz):
    for iat in range(0,len(am)): 
        x[iat] = x[iat] + vx[iat] * dt + ( fx[iat] / (2 * am[iat]) * dt ** 2 )
        y[iat] = y[iat] + vy[iat] * dt + ( fy[iat] / (2 * am[iat]) * dt ** 2 )
        z[iat] = z[iat] + vz[iat] * dt + ( fz[iat] / (2 * am[iat]) * dt ** 2 )
    return(x,y,z)

def update_velocities(dt,am,vx,vy,vz,fx,fy,fz,fx_new,fy_new,fz_new):
    for iat in range(0,len(am)):
        vx[iat] = vx[iat] + ( dt * (fx[iat] + fx_new[iat]) / (2 * am[iat]) )
        vy[iat] = vy[iat] + ( dt * (fy[iat] + fy_new[iat]) / (2 * am[iat]) )
        vz[iat] = vz[iat] + ( dt * (fz[iat] + fz_new[iat]) / (2 * am[iat]) )
     #print(vz[iat])
    return(vx,vy,vz)   

def adjust_velocities(dt,am,vx,vy,vz,fx,fy,fz,fx_new,fy_new,fz_new):
    # at the point where a hop occured, subtract the included forces from an old state and replace them with forces of a new state
    for iat in range(0,len(am)):
        vx[iat] = vx[iat] + ( dt * (-fx[iat] + fx_new[iat]) / (2 * am[iat]) )
        vy[iat] = vy[iat] + ( dt * (-fy[iat] + fy_new[iat]) / (2 * am[iat]) )
        vz[iat] = vz[iat] + ( dt * (-fz[iat] + fz_new[iat]) / (2 * am[iat]) )
     #print(vz[iat])
    return(vx,vy,vz)   
    
def calc_forces(step, at_names, state, nstates, x, y, z, fx, fy, fz,
                pot_eners, ab_initio_file_path):
    """
    Call and collect an external script to calculate ab initio properties (force, energies)
    state = current state - PES for the forcess calc 
    """
    natoms = len(x)

    # Create geom file for which the forces will be calculated
    abinit_geom_file = "abinit_geom.xyz"
    with open (abinit_geom_file, "w") as agf:
         for iat in range(0,len(x)):
             line = ("".join("%2s %2.16e %2.16e %2.16e\n" % (at_names[iat],x[iat]/ANG_BOHR,y[iat]/ANG_BOHR,z[iat]/ANG_BOHR)))
             agf.write(line)
    agf.closed

    # CALL EXTERNAL HANDLING AB-INITIO CALCULATIONS       
    # electronic states starts from 1(ground state) 
    # while code starts from 0 index due to python 
    state = state + 1 
    abinit_inputs = "{} {}  {}  {} {} {}".format(ab_initio_file_path, abinit_geom_file, natoms, state, nstates, step)

    try:
        abinit_proc = subprocess.run(abinit_inputs, stdout=None, stderr=subprocess.PIPE, shell = True, check = True)	
    except subprocess.CalledProcessError as cpe: 
        print("Return code: {}\nError: {}".format(cpe.returncode, cpe.stderr))
        error_exit(4)

    # COLLECT DATA (TODO: diabatization - needs more forces)
    if  re.search(r'g09',ab_initio_file_path):
            grad = 1   # gaus exports forces -1
    elif re.search(r'molpro',ab_initio_file_path):
            grad = -1  #molpro export gradient
    with open ("gradients.dat", "r") as gef:
           
        for st in range(0, nstates):
            pot_eners[st] = float(gef.readline()) 
        for iat in range(0, natoms):
            line = gef.readline().split(" ")
            # FX FY FZ, gradient to forces 
            # TO DO Gaussian has forces, molpro gradients
            fx[iat] = grad * np.float64(line[0])
            fy[iat] = grad * np.float64(line[1])
            fz[iat] = grad * np.float64(line[2])    
    gef.closed
    return(fx , fy, fz, pot_eners)

def calc_energies(
    step, time, natoms, am, state, pot_eners, 
    vx, vy, vz, Etot_init, Etot_prev, ener_thresh):
     
    Ekin = 0.000
    for iat in range(0,natoms):
        vv = vx[iat] ** 2 + vy[iat] ** 2 + vz[iat] ** 2
        Ekin = Ekin + (0.5 * am[iat] * vv)
        Epot = pot_eners[state]  # state 0(GS), 1 (1.ex. state),..
        Etot = Ekin + Epot
        dE = (Etot - Etot_init)
        if (dE * AU_EV) >= ener_thresh: 
            print("Total energy change since start {} ".format(dE * AU_EV),
                  "larger then threshold {}.".format(ener_thresh))
            error_exit(7)
        dE_step = (Etot - Etot_prev)
    print_energies(step, time, Ekin, Epot, Etot, dE, dE_step)
    print_pes(time, step, pot_eners)

    return(Ekin, Epot, Etot, dE, dE_step)

def rescale_velocities(vx, vy, vz, v_scaling_fac):

    vx = [xx * v_scaling_fac for xx in vx] 
    vy = [yy * v_scaling_fac for yy in vy] 
    vz = [zz * v_scaling_fac for zz in vz] 

    return(vx, vy, vz)
    # Windows installed ubuntu has rather complicated path
    # if re.search(r'win',sys.platform):
    # testpath = "/mnt/c/Users/chalabaj/Documents/Coding/snafu-master/ABINITIO/test.sh" 
    # abinit_inputs = "wsl {} {}  {}  {}  {}".format(testpath, abinit_geom_file, natoms, state, nstates, step)
    # elif re.search(r'linux',sys.platform):  
