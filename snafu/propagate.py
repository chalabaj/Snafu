import math
import sys, os
import numpy as np
import random
import time

def verlet_step(natoms,dt,am,x,y,z,vx,vy,vz):
      

    for iat in range(0,natoms):
      
   # update positions: Rn(t + Δt) = Rn(t) + Δt*Vn(t) +Δt^2/(2Mn)*Fn(t)  - for n atom at time t
      x[iat] = x[iat] + vx[iat] * dt + 1/(2*am[iat]) * fx[iat] * (dt**2)
      y[iat] = y[iat] + vy[iat] * dt + 1/(2*am[iat]) * fy[iat] * (dt**2)
      z[iat] = z[iat] + vz[iat] * dt + 1/(2*am[iat]) * fz[iat] * (dt**2)
      print(iat,"x:",x[iat],"fx:",fx[iat]/am[iat])
   
    # calc_force: Fn(t+dt) F_new
     
   # update_velocities: Vn(t+dt) Vn(t + Δt) = Vn(t) + Δt/(2Mn)*(Fn(t) + Fn(t + Δt))
    for iat in range(0,natoms):
      
      # update_velocities: Vn(t + Δt) = Vn(t) + Δt/(2Mn)*(Fn(t) + Fn(t + Δt))
      vx[iat] = vx[iat] + 1/(2*am[iat]) * (fx[iat] + fx_new[iat])
      vy[iat] = vy[iat] + 1/(2*am[iat]) * (fy[iat] + fy_new[iat])
      vz[iat] = vz[iat] + 1/(2*am[iat]) * (fz[iat] + fz_new[iat])
    
    return(x,y,z,fx_new,fy_new,fz_new, pot_eners)   # F are needed for diabatization scheme, potential energizi for APLZ scheme
