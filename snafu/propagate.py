import math
import sys, os
import numpy as np
import random
import time

def verlet_step(natoms,dt,am,x,vx,fx):
    
    for iat in range(0,natoms):
      x[iat] = x[iat] + vx[iat]*dt + fx[iat]/am[iat]*(dt**2)/2
      print(iat,"x:",x[iat],"fx:",fx[iat]/am[iat])
    return x    
    
    #x_new = x + vx*dt + a*(dt**2)/2
	 #v_new = v + (a(r) + a(r_new))/2 * dt
    #print(x_new)
    #return(x_new)#y_new,z_new,vx_new,vy_new,vz_new,e_new)    

#return x_new, y_new, z_new, v_new, v_new, v_new,
    
"""
        
    Main rutine to propagate atoms
    Calling external code to calc. force, energies
    Calcul transition probabilities 
    Readjust velocities
    Check. energy conservation
"""