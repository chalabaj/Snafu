
ang_bohr = 1.889726132873e0     # agstroms to bohrs

def print_positions(step,time,natoms, at_names, x, y, z):
    bohr_ang =1/ang_bohr     # bohr to ang units
    
    with open ("movie.xyz", "a") as mov:
     header = ("{} \n".format(natoms))
     mov.write(header)
     comment = ("Step: {}      Time_fs:{}\n".format(step,time))
     mov.write(comment)
     for iat in range(0,natoms):     
      line = ("".join("%2s %5.8f %5.8f %5.8f\n"  %(at_names[iat],x[iat]*bohr_ang,y[iat]*bohr_ang,z[iat]*bohr_ang)))
      mov.write(line)
    mov.closed
    return()

def print_velocities(step,time,natoms, at_names, vx, vy, vz):
     
    with open ("velocities.xyz", "a") as vel:
     header = ("{} \n".format(natoms))
     vel.write(header)
     comment = ("Step: {}      Time_fs:{}\n".format(step,time))
     vel.write(comment)
     for iat in range(0,natoms):
       line = ("".join("%2s %5.8f %5.8f %5.8f\n"  %(at_names[iat],vx[iat],vy[iat],vz[iat])))
       vel.write(line)
    vel.closed
    return()
 
def print_energies(step,time,Ekin,Epot,Etot,dE):
    with open ("energies.dat", "a") as ef:
        if step == 0:
            headline = "# Time,  Ekinetic/au,  Epotential/au,  Etotal/au,  dE/eV\n"
            ef.write(str(headline))
        line = "{:10.10f} {:20.10f} {:20.10f} {:20.10f} {:20.10f}\n".format(time,Ekin,Epot,Etot,dE)
        ef.write(str(line))
    ef.closed
    return()  
    
def print_pes(time, step, pot_eners):
    with open ("PES.dat", "a") as pesf:
        if step == 0:
            headline = "# Time,  E(GS)/au,  E(1. ex)/au,....\n"
            pesf.write(str(headline))

        line = ("{:10.10f}".format(time)
                + ' '.join('{:20.10f}'.format(pot_eners[st]) 
                for st in range(0, len(pot_eners)))
                + "\n")
        pesf.write(str(line))
    pesf.closed
    return()
    
def print_snafu():
   print("   SSSSSSSSSSSSSSS NNNNNNNN        NNNNNNNN               AAA               FFFFFFFFF FFFFFFFFFFFUUUUUUUU     UUUUUUUU")
   print( " SS:::::::::::::::SN:::::::N       N::::::N              A:::A              F:::::::::::::::::::FU::::::U     U::::::U")
   print( "S:::::SSSSSS::::::SN::::::::N      N::::::N             A:::::A             F:::::::::::::::::::FU::::::U     U::::::U")
   print( "S:::::S     SSSSSSSN:::::::::N     N::::::N            A:::::::A            FF:::::FFFFFFFFF::::FUU:::::U     U:::::UU")
   print( "S:::::S            N::::::::::N    N::::::N           A:::::::::A             F::::F       FFFFFF U:::::U     U:::::U ")
   print( "S:::::S            N:::::::::::N   N::::::N          A:::::A:::::A            F::::F              U:::::D     D:::::U ")
   print( " S::::SSSS         N:::::::N::::N  N::::::N         A:::::A A:::::A           F:::::FFFFFFFFFF    U:::::D     D:::::U ")
   print( "  SS::::::SSSSS    N::::::N N::::N N::::::N        A:::::A   A:::::A          F::::::::::::::F    U:::::D     D:::::U ")
   print( "    SSS::::::::SS  N::::::N  N::::N:::::::N       A:::::A     A:::::A         F::::::::::::::F    U:::::D     D:::::U ")
   print( "       SSSSSS::::S N::::::N   N:::::::::::N      A:::::AAAAAAAAA:::::A        F:::::FFFFFFFFFF    U:::::D     D:::::U ")
   print( "            S:::::SN::::::N    N::::::::::N     A:::::::::::::::::::::A       F::::F              U:::::D     D:::::U ")
   print( "            S:::::SN::::::N     N:::::::::N    A:::::AAAAAAAAAAAAA:::::A      F::::F              U::::::U   U::::::U ")
   print( "SSSSSSS     S:::::SN::::::N      N::::::::N   A:::::A             A:::::A   FF::::::FF            U:::::::UUU:::::::U ")
   print( "S::::::SSSSSS:::::SN::::::N       N:::::::N  A:::::A               A:::::A  F:::::::FF             UU:::::::::::::UU  ")
   print( "S:::::::::::::::SS N::::::N        N::::::N A:::::A                 A:::::A F:::::::FF               UU:::::::::UU    ")
   print( " SSSSSSSSSSSSSSS   NNNNNNNN         NNNNNNNAAAAAAA                   AAAAAAAFFFFFFFFFF                 UUUUUUUUU      \n")
   return()
