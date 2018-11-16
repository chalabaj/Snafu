import os
import numpy as np
import shutil
try:
    from constants import *
    from errors import error_exit
except ImportError as ime:
    # module could have been removed or module file renamed
    if ime.name is None:  
        print("Import in some of the modules ({})".format(ime.name),
              "in snafu dir failed. Exiting...")
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

def print_positions(step,sim_time,natoms, at_names, x, y, z, restart):
    
    with open ("movie.xyz", "a") as mov:
     header = ("{} \n".format(natoms))
     mov.write(header)
     comment = ("Step: {}      Time_fs: {}\n".format(step,sim_time))
     mov.write(comment)
     for iat in range(0, natoms):     
      line = ("".join("%2s %5.8f %5.8f %5.8f\n"  %(at_names[iat],x[iat]*BOHR_ANG,y[iat]*BOHR_ANG,z[iat]*BOHR_ANG)))
      mov.write(line)
    mov.closed
    return()

def print_velocities(step,sim_time,natoms, at_names, vx, vy, vz, restart):
    with open ("velocities.xyz", "a") as vel:
     header = ("{} \n".format(natoms))
     vel.write(header)
     comment = ("Step: {}      Time_fs:{}\n".format(step,sim_time))
     vel.write(comment)
     for iat in range(0, natoms):
       line = ("".join("%2s %5.8f %5.8f %5.8f\n"  %(at_names[iat],vx[iat],vy[iat],vz[iat])))
       vel.write(line)
    vel.closed
    return()

def print_energies(step, sim_time, Ekin, Epot, Etot, dE, dE_step, restart):
    with open ("energies.dat", "a") as ef:
        if step == 1:
            headline = "# Time,  Ekin,  Epot,  Etot,  dE,  dE_step (a.u.)\n"
            ef.write(str(headline))
            dE = 0.0
            dE_step  = 0.0
        line = "{:>10.2f} {:20.10f} {:20.10f} {:20.10f} {:20.10f}  {:20.10f}\n".format(sim_time,Ekin,Epot,Etot,dE, dE_step )
        ef.write(str(line))
    ef.closed
    return()  

def print_pes(sim_time, step, pot_eners):
    with open ("PES.dat", "a") as pesf:
        if step == 1:
            headline = "# Time,  E(GS)/au,  E(1. ex)/au,....\n"
            pesf.write(str(headline))

        line = ("{:10.10f}".format(sim_time)
                + ' '.join('{:20.10f}'.format(pot_eners[st]) 
                for st in range(0, len(pot_eners)))
                + "\n")
        pesf.write(str(line))
    pesf.closed
    return()

def print_state(step, sim_time, state):
    with open ("state.dat", "a") as stf:
        if step == 1:
            headline = "# Time,  Electronic state(0 = gs)\n"
            stf.write(str(headline)) 
        line = ("{:7.4f} {:4d}\n".format(sim_time, state)) 
        stf.write(str(line))   
    stf.closed
    return()

def print_snafu():
   print("   SSSSSSSSSSSSSSS NNNNNNNN        NNNNNNNN                AAA               FFFFFFFFF FFFFFFFFFFFUUUUUUUU     UUUUUUUU")
   print( " SS:::::::::::::::SN:::::::N       N::::::N              A:::A              F:::::::::::::::::::FU::::::U     U::::::U")
   print( "S:::::SSSSSS::::::SN::::::::N      N::::::N             A:::::A             F:::::::::::::::::::FU::::::U     U::::::U")
   print( "S:::::S     SSSSSSSN:::::::::N     N::::::N            A:::::::A            FF:::::FFFFFFFFF::::FUU:::::U     U:::::UU")
   print( "S:::::S            N::::::::::N    N::::::N           A:::::::::A             F::::F       FFFFFF U:::::U     U:::::U ")
   print( "S:::::S            N:::::::::::N   N::::::N          A:::::A:::::A            F::::F              U:::::D     D:::::U ")
   print( " S::::SSSS         N:::::::N::::N  N::::::N         A:::::A A:::::A           F:::::FFFFFFFFFF    U:::::D     D:::::U ")
   print( "  SS::::::SSSSS    N::::::N N::::N N::::::N        A:::::A   A:::::A          F::::::::::::::F    U:::::D     D:::::U ")
   print( "    SSS::::::::SS  N::::::N  N::::N:::::::N       A:::::A *** A:::::A         F::::::::::::::F    U:::::D     D:::::U ")
   print( "       SSSSSS::::S N::::::N   N:::::::::::N      A:::::AAAAAAAAA:::::A        F:::::FFFFFFFFFF    U:::::D     D:::::U ")
   print( "            S:::::SN::::::N    N::::::::::N     A:::::::::::::::::::::A       F::::F              U:::::D     D:::::U ")
   print( "            S:::::SN::::::N     N:::::::::N    A:::::AAAAAAAAAAAAA:::::A      F::::F              U::::::U   U::::::U ")
   print( "SSSSSSS     S:::::SN::::::N      N::::::::N   A:::::A             A:::::A   FF::::::FF            U:::::::UUU:::::::U ")
   print( "S::::::SSSSSS:::::SN::::::N       N:::::::N  A:::::A               A:::::A  F:::::::FF             UU:::::::::::::UU  ")
   print( "S:::::::::::::::SS N::::::N        N::::::N A:::::A                 A:::::A F:::::::FF               UU:::::::::UU    ")
   print( " SSSSSSSSSSSSSSS   NNNNNNNN         NNNNNNNAAAAAAA                   AAAAAAAFFFFFFFFFF                 UUUUUUUUU      \n")
   return()
