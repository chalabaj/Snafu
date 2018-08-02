try:
    import sys
    import os
    SNAFU_EXE = os.environ['SNAFU_DIR']
    sys.path.append(os.path.join(SNAFU_EXE, "snafu"))
    sys.path.append(SNAFU_EXE)
    from landauzener import (
        calc_lz_hopp
    )
    import numpy as np
    from constants import *
except ImportError as ime:
    if ime.name is None:  # module could have been removed or module file renamed
        print("Import in some of the modules ({}) in snafu dir failed. Exiting...".format(ime.name))
        exit(1)
    else:
        print("Module {} not found.".format(ime.name),
              "Make sure that {} contains snafu folder".format(SNAFU_EXE),
              "with: {}".format('\n'.join(modules_files)))
        exit(1)
except KeyError as ke:
    print("SNAFU_DIR is not set.",
          "See '.bashrc' file in your home directory",
          "or use 'env' command and make sure $SNAFU_DIR is exported.",
          "\nHint: export SNAFU_DIR=/path/to/SNAFU")
    exit(1)
else:
    print("All modules loaded succesfully. Starting...\n \n")

step = 0
nstates = 8
init_state = 2
Ekin = 0.1   #ua
dt = 1
pot_eners = np.zeros(nstates,dtype=np.float64) 
#print(pot_eners)
step = 0
nhops = 0
inputfile = sys.argv[1]
outfile = "{}_hopping-initstate_{}-{}auEkin_dt{}.dat".format(inputfile,init_state,Ekin,dt)
os.system("rm {}".format(outfile))
baseline = 196.9
with open(inputfile,'r') as gf:  
      
     gf.readline()  # header
     state = init_state
     for l in gf:
         step += 1
         print(step)
         line = gf.readline().split(" ")
         for st in range(0, nstates):
             pot_eners[st] = np.float64(line[st])   
         #print(pot_eners)
         if step == 1:
             pot_eners_array = np.copy(pot_eners)
             #print(pot_eners_array)
         elif step == 2:
             pot_eners_array = np.vstack((pot_eners_array, pot_eners))
             #print(pot_eners_array)
         else:
             hop, outstate, pot_eners_array, prob = calc_lz_hopp("l",
                                                                 state ,
                                                                 pot_eners,
                                                                 pot_eners_array,
                                                                 Ekin, dt)     
             
             if  hop:
                 state = outstate
                 nhops += 1
             
             with open(outfile,'a') as of: 
                 #line = ("{} {}".format(step, state) )

                 line = ("{}  {} {} ".format(step, 0.1*state-baseline, prob-baseline)
                         + ' '.join('{:5.5f}'.format(pot_eners_array[1][st]) 
                         for st in range(0, nstates))
                         + "\n")
                 of.write(line)
             of.close()
     print("Number of hops: {}".format(nhops))
     
gf.close()
