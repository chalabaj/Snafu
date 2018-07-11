# -*- coding: utf-8 -*-
"""
Created on Tue Jul 10 14:07:22 2018

@author: chalabaj
"""


def loop_step():
    """
    Main rutine to propagate atoms
    Calling external code to calc. force, energies
    Calcul transition probabilities 
    Readjust velocities
    Check. energy conservation
    """
#  x_old  x(t-dt)
#  x_     x(t)
#  x_new  x(t+dt)
    #cacl_force
    for step in range (o,maxsteps):
        print("verlet")

        
             
        
    return

def verlet_step(x,y,z,m,t):
    for iat in range(1,natoms+1):   # upper index excluded
     print("for atom")   
    #eq 1 px(iat,t + 1/2 dt) = px(iat,t) + 1/2*dt*fi(iatt))
    #eq 2
    #call forces
    # eq 3
    return()


#return x_new, y_new, z_new, v_new, v_new, v_new,