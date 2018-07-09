"""
SNAFU code
SNAFU can do ab initio nonadiabatic MD withot calculating non-adiabatic couplings
Two scheme with/without diabatization

1) adiabatic-potential-based formula derived within Landau-Zener:
AK Belyaev, PHYSICAL REVIEW A 84, 014701 (2011),DOI: 10.1103/PhysRevA.84.014701
AK Belyaev, The Journal of Chemical Physics 147, 234301 (2017); doi: 10.1063/1.5000718

2) diabatization scheme: Le Yu, Phys.Chem.Chem.Phys., 2014, 16, 25883
"""

import sys
import numpy as np

# Constants
au_fs = 0.02
au_eV = 27.21139

# Observed variable
aEk = 0.0   # average kinetic energy
cPE = 0.0   # current potential energy to save
cKE = 0.0   # current kinetic energy to save

#  x_old  x(t-dt)
#  x_     x(t)
#  x_new  x(t+dt)

if __name__ == "__init__":
    print_info(start)
    # file_check
    # read input
    # create array  - pos, velocities, energies, gradients
    # load initial pos, velocities

def main_loop():
    
    #cacl_force
    for step in range (o,maxsteps):
        #verlet_step
        #calc_hop
        # vel_adjustment
        #calc_energies
        #print_info(step, pos, ener)
        
             
        
    return

def verlet_step(x,y,z,m,t):
    for iat in range(1,natoms+1)   # upper index excluded
        
    #eq 1 px(iat,t + 1/2 dt) = px(iat,t) + 1/2*dt*fi(iatt))
    #eq 2
    #call forces
    # eq 3
return x_new, y_new, z_new, v_new, v_new, v_new,