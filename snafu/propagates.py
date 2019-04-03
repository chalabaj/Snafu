"""
This module contain MD integrator parts: update velocities and position.
Function "calc_forces" calls external ab initio program to calculate forces and then read the forces
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
current_module = sys.modules[__name__]
try:
    from errors import error_exit
    from constants import *
    from tera_propagates import (
        receive_tera, send_tera, exit_tera
    )
except ImportError as ime:
    print("Module {} in {} not found.".format(ime,current_module))
    exit(1)

def update_positions(
        dt, am, x, y, z, x_new, y_new, z_new, vx, vy, vz, fx, fy, fz):
    for iat in range(0,len(am)):
        x_new[iat] = x[iat] + vx[iat]*dt + ( fx[iat]/(2*am[iat])*dt**2 )
        y_new[iat] = y[iat] + vy[iat]*dt + ( fy[iat]/(2*am[iat])*dt**2 )
        z_new[iat] = z[iat] + vz[iat]*dt + ( fz[iat]/(2*am[iat])*dt**2 )
    return(x_new, y_new, z_new)

def update_velocities(dt, am, vx, vy, vz, fx, fy, fz, fx_new, fy_new, fz_new):
    vxx = np.copy(vx)
    vyy = np.copy(vy)
    vzz = np.copy(vz)
    for iat in range(0,len(am)):
        vxx[iat] = vx[iat] + ( dt*(fx[iat] + fx_new[iat])/(2*am[iat]) )
        vyy[iat] = vy[iat] + ( dt*(fy[iat] + fy_new[iat])/(2*am[iat]) )
        vzz[iat] = vz[iat] + ( dt*(fz[iat] + fz_new[iat])/(2*am[iat]) )
    return(vxx, vyy, vzz)   

def adjust_velocities(dt, am, vx, vy, vz, fx, fy, fz, fx_new, fy_new, fz_new):
    vxx = np.copy(vx)
    vyy = np.copy(vy)
    vzz = np.copy(vz)
    # at the point where a hop occured, 
    # subtract the included forces from an old state 
    # and replace them with forces of a new state
    for iat in range(0,len(am)):
        vxx[iat] = vx[iat] + ( dt*(-fx[iat] + fx_new[iat])/(2*am[iat]) )
        vyy[iat] = vy[iat] + ( dt*(-fy[iat] + fy_new[iat])/(2*am[iat]) )
        vzz[iat] = vz[iat] + ( dt*(-fz[iat] + fz_new[iat])/(2*am[iat]) )
    return(vxx, vyy, vzz)   
    
def calc_forces(step, at_names, state, nstates, x, y, z, fx_new, fy_new, fz_new, 
                pot_eners, ab_initio_file_path, 
                tera_mpi, comm, sim_time, MO, CiVecs, NAC, blob, SMatrix, civec_size, nbf_size, blob_size, qmcharges, TDip, Dip):
    """ 
    Call and collect an external script to calculate ab initio properties (force, energies)
    state = current state - PES for the forces to calc 
    """
    natoms = len(x)
    grad = -1
    if not tera_mpi:
        # GRADIENT/FORCES check
        if re.search(r'g09', ab_initio_file_path) or re.search(r'gaus09', ab_initio_file_path):
                grad = 1   # gaus exports forces -1
        elif re.search(r'molpro', ab_initio_file_path):
                grad = -1  # molpro exports gradients
                state = state + 1 #molpro index starts from 1
        elif re.search(r'orca', ab_initio_file_path):
                grad = -1  #orca exports gradients
        elif re.search(r'tera', ab_initio_file_path):
                grad = -1  #tera exports gradients                    
        # Create geom file for which the forces will be calculated
        abinit_geom_file = "abinit_geom.xyz"
        with open (abinit_geom_file, "w") as agf:
             for iat in range(0,natoms):
                 line = ("".join("%2s %2.16e %2.16e %2.16e\n" % (at_names[iat],x[iat]/ANG_BOHR,y[iat]/ANG_BOHR,z[iat]/ANG_BOHR)))
                 agf.write(line)
        agf.closed
    
        # CALL EXTERNAL AB-INITIO CALCULATIONS       
        abinit_inputs = "{} {}  {}  {} {} {}".format(ab_initio_file_path, abinit_geom_file, natoms, state, nstates, step)
    
        try:
            abinit_proc = subprocess.run(abinit_inputs, stdout=None, stderr=subprocess.PIPE, shell = True, check = True)	
        except subprocess.CalledProcessError as cpe: 
            print("Return code: {}\nError: {}".format(cpe.returncode, cpe.stderr))
            error_exit(4)
    
        with open ("gradients.dat", "r") as gef:
            for st in range(0, nstates):
                pot_eners[st] = float(gef.readline()) 
            for iat in range(0, natoms):
                line = gef.readline().split(" ")
                # FX FY FZ, gradient to forces 
                fx_new[iat] = grad*np.float64(line[0])
                fy_new[iat] = grad*np.float64(line[1])
                fz_new[iat] = grad*np.float64(line[2])    
        gef.closed
    else:
        try:
            send_tera(comm, natoms, nstates, state, sim_time, x ,y, z, MO, CiVecs, blob, civec_size, nbf_size, blob_size)
            
            fx_new, fy_new, fz_new, pot_eners, MO, CiVecs, blob = receive_tera(comm, natoms, nstates, state, pot_eners, fx_new, fy_new, fz_new,
                                                                               MO, CiVecs, blob, SMatrix, NAC, TDip, Dip, qmcharges, civec_size, nbf_size, blob_size)

            fx_new = grad*fx_new
            fy_new = grad*fy_new
            fz_new = grad*fz_new
        except Exception as excpt:
            print("Something went wrong during MPI SEND/RECEIVE.",
                  "\n{}".format(excpt))
            finish_tera(comm)
            error_exit(error_exit(15, str("Error during sending/receive TC data {}".format(excpt)))) 
        
    return(fx_new , fy_new, fz_new, pot_eners, MO, CiVecs, blob)

def calc_energies(
    step, time, natoms, am, state, pot_eners, 
    vx, vy, vz, Etot_init, Etot_prev, ener_thresh):

    Ekin = 0.000

    for iat in range(0,natoms):
        vv = vx[iat]**2 + vy[iat]**2 + vz[iat]**2
        Ekin = Ekin + (0.5*am[iat]*vv)

    Epot = pot_eners[state]  # state 0(GS), 1 (1.ex. state),..
    Etot = Ekin + Epot
    dE = (Etot - Etot_init)
    dE_step = (Etot - Etot_prev)

    if abs(dE * AU_EV) >= ener_thresh and step > 1: 
        print("Total energy change since start {} ".format(dE * AU_EV),
              "larger then threshold {}.".format(ener_thresh))
        error_exit(7)

    return(Ekin, Epot, Etot, dE, dE_step)

def rescale_velocities(vx, vy, vz, v_scaling_fac):
    return(vx*v_scaling_fac, vy*v_scaling_fac, vz*v_scaling_fac)
