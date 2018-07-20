import math
import sys, os
import subprocess
import random
import time
from snafu.errors import error_exit
import re
# CONSTANTS
au_fs = 0.02418884326505      #atomic units to femtosecs
au_eV = 27.21139
amu   = 1822.888484264545             # atomic mass unit  me = 1 AMU*atomic weight
ang_bohr = 1.889726132873     # agstroms to bohrs
bohr_ang = 1/ang_bohr     # bohr to ang units
#====================================================================================================


def update_positions(natoms,dt,am,x,y,z,vx,vy,vz,fx,fy,fz):
# update positions: Rn(t + Δt) = Rn(t) + Δt*Vn(t) +Δt^2/(2Mn)*Fn(t)  - for n atom at time t
    #print(x,y,z)    print(vx,vy,vz)    print(fx,fy,fz)

    for iat in range(0,natoms):   
      x[iat] = x[iat] + vx[iat] * dt + ( 0.5 * fx[iat] * (dt ** 2) / (am[iat]) )
      y[iat] = y[iat] + vy[iat] * dt + ( 0.5 * fy[iat] * (dt ** 2) / (am[iat]) )
      z[iat] = z[iat] + vz[iat] * dt + ( 0.5 * fz[iat] * (dt ** 2) / (am[iat]) )
      #print(x[iat],vx[iat] * dt,(1/(2*am[iat]) * fx[iat] * (dt**2)))
    return(x,y,z)


def update_velocities(natoms,dt,am,vx,vy,vz,fx,fy,fz,fx_new,fy_new,fz_new):
# update_velocities: Vn(t+dt) Vn(t + Δt) = Vn(t) + Δt/(2Mn)*(Fn(t) + Fn(t + Δt))

    for iat in range(0,natoms):
     vx[iat] = vx[iat] + ( 0.5 * dt * (fx[iat] + fx_new[iat]) / am[iat] )
     vy[iat] = vy[iat] + ( 0.5 * dt * (fy[iat] + fy_new[iat]) / am[iat] )
     vz[iat] = vz[iat] + ( 0.5 * dt * (fz[iat] + fz_new[iat]) / am[iat] )
     #print(vz[iat])
    return(vx,vy,vz)   
    

def calc_forces(natoms, at_names, state, nstates, ab_initio_file_path, x, y, z, fx_new, fy_new, fz_new, pot_eners):
    """
    Call and collect an external script to calculate ab initio properties (force, energies)
    state = current state - PES for the forcess calc 
    """
    # Create geom file for which the forces will be calculated
    abinit_geom_file = "abinit_geom.xyz"
    with open (abinit_geom_file, "w") as agf: #ab init geom file
         for iat in range(0,natoms):
             line = ("".join("%2s %4.18e %4.18e %4.18e\n"  %(at_names[iat],x[iat]*bohr_ang,y[iat]*bohr_ang,z[iat]*bohr_ang)))
             agf.write(line)
    agf.closed
    
    # Windows installed ubuntu has rather complicated path
    testpath = "/mnt/c/Users/chalabaj/Documents/Coding/snafu-master/ABINITIO/test.sh" 
    
    state = state + 1 # ab initio code starts from 1(ground state) while code starts from 0 index due to python 
    
    if re.search(r'win',sys.platform):   
       abinit_inputs = "wsl {} {}  {}  {}  {}".format(testpath, abinit_geom_file, natoms, state, nstates)
    elif re.search(r'linux',sys.platform): 
       abinit_inputs = "{} {}  {}  {}  {}".format(ab_initio_file_path, abinit_geom_file, natoms, state, nstates)
       
    # CALL EXTERNAL SCRIPT WHICH WILL HANDLE AB INITIO CALCULATIONS       
    try:
        abinit_proc = subprocess.run(abinit_inputs, stdout= None, stderr= subprocess.PIPE, shell = True, check = True)	
    except subprocess.CalledProcessError as cpe: 
        print("Return code: {}\nError: {}".format(cpe.returncode, cpe.stderr))
        error_exit(4)
    #else:
     #   print("Calculating forces.")    # {}".format(abinit_proc.stdout))
        
    # Check return status of the process:     # not necesarry  exception handle this: if abinit_call.returncode = 0:     error_exit(4)   
    
    # COLLECT DATA
    #TODO: diabatization - need more forces
    with open ("gradients.dat", "r") as gef:   # gradient energy file
     for st in range(0,nstates):
        pot_eners[st] = float(gef.readline())  # comment 
     for iat in range(0,natoms):
        line = gef.readline().split()
        fx_new[iat] = float(line[0])
        fy_new[iat] = float(line[1])
        fz_new[iat] = float(line[2])   #  X Y Z format for each atoms
       
    gef.closed
    return(fx_new , fy_new, fz_new, pot_eners)

def calc_energies(step, time, natoms, am, state, pot_eners, vx, vy, vz, Etot_init):
    Ekin = 0.0
    for iat in range(0,natoms):

         Ekin = Ekin + (0.5*am[iat] * (vx[iat] **2 + vy[iat]**2 + vz[iat]**2))
         #print(Ekin)
         Epot = pot_eners[state] #state 0(GS), 1 (1.ex. state),.....
         Etot = Ekin + Epot
         dE = (Etot - Etot_init) 

    with open ("energies.dat", "a") as ef:
       if step == 0:
        line = "# Time,  Ekinetic/au,  Epotential/au,  Etotal/au,  dE/eV\n"
        ef.write(str(line))
        #line = "{:10.10f} {:20.10f} {:20.10f} {:20.10f} {:20.10f}\n".format(time,Ekin,Epot,Etot,dE*au_eV)
        #ef.write(str(line))
       else:
        line = "{:10.10f} {:20.10f} {:20.10f} {:20.10f} {:20.10f}\n".format(time,Ekin,Epot,Etot,dE*au_eV)
       ef.write(str(line))
    ef.closed

    return(Ekin,Epot,Etot,dE)

    