import numpy as np

def read_restart():
    #with open restart.in as rsf:
    np.loadtxt("restart.in", dtype=np.float64, delimiter=None, skiprows=8)  
    return()
    
if __name__ == "__main__":
    read_restart()
    
    exit(0)
    
