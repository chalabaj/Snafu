import numpy as np

def read_restart():
    #with open restart.in as rsf:
    rst_file = "restart.in"
    with open(rst_file, 'r') as rstf:
        step = rstf.readline().split()[1]
        print(step)
    rstf.closed
    a = np.loadtxt(rst_file, dtype=np.float64, delimiter=None, skiprows=8)  
    print(a)
    return()
    
if __name__ == "__main__":
    read_restart()
    
    exit(0)
    
