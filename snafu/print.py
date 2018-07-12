import numpy as np

def print_positions(step,time,natoms, at_names, x, y, z):
    
    with open ("movie.xyz", "w") as mov:
     header = ("{} \n".format(natoms))
     mov.write(header)
     comment = ("Step: {}      Time_fs:{}".format(step,time))
     mov.write(comment)
     for iat in range(0,natoms):
      
      line = ("".join("%2s %2.4f %2.4f %2.4f"  %(at_names[iat],x[iat],y[iat],z[iat])))
      print(line)
      mov.write(line)
    mov.closed
    return()

