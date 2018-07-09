"""
main rutine
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
