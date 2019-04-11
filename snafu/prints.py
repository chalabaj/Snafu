import os
import sys
try:
    from constants import *
    from errors import error_exit   
except ImportError as ime:
    error_exit(19, "Module {} in {} not found.".format(ime,current_module))
try:
    tera_mpi = int(os.environ['MPI_TERA'])
except KeyError as ke:
    print("MPI_TERA variable was not exported, assuming MPI_TERA=0. Warning: this may cause deadlock if MPI has been already initiated")
    tera_mpi = 0
if tera_mpi:
    from tera_propagates import global_except_hook
    sys.excepthook = global_except_hook  

# ---------------------------------------------------------------------------------------------------------------------------------------------    

def print_positions(step,sim_time,natoms, at_names, x, y, z, mov_file):
    mov = mov_file
    header = ("{} \n".format(natoms))
    mov.write(header)
    comment = ("Step: {}      Time_fs: {}\n".format(step, sim_time))
    mov.write(comment)
    xx = (x*BOHR_ANG).tolist()
    yy = (y*BOHR_ANG).tolist()
    zz = (z*BOHR_ANG).tolist()
    for iat in range(0, natoms):
        line = "{} {:12.8f} {:12.8f} {:12.8f}\n".format(at_names[iat], xx[iat],
                                                        yy[iat], zz[iat])
        mov.write(line)
    mov.flush()
    return()

def print_velocities(step,sim_time,natoms, at_names, vx, vy, vz, vel_file):
    vel = vel_file
    header = ("{} \n".format(natoms))
    vel.write(header)
    comment = ("Step: {}      Time_fs:{}\n".format(step, sim_time))
    vel.write(comment)
    for iat in range(0, natoms):
       line = ("".join("%2s %5.8f %5.8f %5.8f\n"  %(at_names[iat], vx[iat], vy[iat], vz[iat])))
       vel.write(line)
    vel.flush()
    return()

def print_energies(step, write_freq, sim_time, Ekin, Epot, Etot, dE, dE_step, eners_file):
    ef = eners_file
    if step == 1 or step == write_freq:
        headline = "# Time_fs,  Ekin,  Epot,  Etot,  dE,  dE_step (a.u.)\n"
        ef.write(str(headline))
    line = "{:>10.4f} {:20.10f} {:20.10f} {:20.10f} {:20.10f}  {:20.10f}\n".format(sim_time,Ekin,Epot,Etot,dE, dE_step )
    ef.write(str(line))
    ef.flush()
    return()  

def print_pes(step, write_freq, sim_time, pot_eners, pes_file):
    if step == 1 or  step == write_freq:
        headline = "# Time_fs,  E(GS)/au,  E(1. ex)/au,....\n"
        pes_file.write(str(headline))
    line = ("{:10.4f}".format(sim_time)
                + ' '.join('{:20.10f}'.format(pot_eners[st]) 
                for st in range(0, len(pot_eners)))
                + "\n")
    pes_file.write(str(line))
    pes_file.flush()
    return()

def print_state(step, write_freq, sim_time, state, state_file):
    if step == 1 or  step == write_freq:
        headline = "# Time_fs,  Electronic state(0 = gs)\n"
        state_file.write(str(headline)) 
    line = ("{:10.4f} {:4d}\n".format(sim_time, state)) 
    state_file.write(str(line))   
    state_file.flush()
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
