
def print_positions(step,time,natoms, at_names, x, y, z):
    
    #line = "  ".join("%1s %1d" "%1d" "%1d" %(x[1], y[1], z[1]))
    with open ("movie.xyz", "w") as mov:
    
     mov.write(natoms,"\n")
     mov.write("Step: ",step,"     ","Time_fs: ",time,"\n")
     for iat in range(0,natoms):
      
      line = (' {0: }{0:2f.3} {0:2f.3} {0:2f.3}\n'.format(x[iat], y[iat], z[iat]))
      line = "  ".join("%1s %1d" "%1d" "%1d" %(x[iat], y[iat], z[iat]))
     mov.write(line)
    mov.closed
    return()

