
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
    
def print_snafu():
 print("   SSSSSSSSSSSSSSS NNNNNNNN        NNNNNNNN               AAA               FFFFFFFFFFFFFFFFFFFFFFUUUUUUUU     UUUUUUUU")
 print( " SS:::::::::::::::SN:::::::N       N::::::N              A:::A              F::::::::::::::::::::FU::::::U     U::::::U")
 print( "S:::::SSSSSS::::::SN::::::::N      N::::::N             A:::::A             F::::::::::::::::::::FU::::::U     U::::::U")
 print( "S:::::S     SSSSSSSN:::::::::N     N::::::N            A:::::::A            FF::::::FFFFFFFFF::::FUU:::::U     U:::::UU")
 print( "S:::::S            N::::::::::N    N::::::N           A:::::::::A             F:::::F       FFFFFF U:::::U     U:::::U ")
 print( "S:::::S            N:::::::::::N   N::::::N          A:::::A:::::A            F:::::F              U:::::D     D:::::U ")
 print( " S::::SSSS         N:::::::N::::N  N::::::N         A:::::A A:::::A           F::::::FFFFFFFFFF    U:::::D     D:::::U ")
 print( "  SS::::::SSSSS    N::::::N N::::N N::::::N        A:::::A   A:::::A          F:::::::::::::::F    U:::::D     D:::::U ")
 print( "    SSS::::::::SS  N::::::N  N::::N:::::::N       A:::::A     A:::::A         F:::::::::::::::F    U:::::D     D:::::U ")
 print( "       SSSSSS::::S N::::::N   N:::::::::::N      A:::::AAAAAAAAA:::::A        F::::::FFFFFFFFFF    U:::::D     D:::::U ")
 print( "            S:::::SN::::::N    N::::::::::N     A:::::::::::::::::::::A       F:::::F              U:::::D     D:::::U ")
 print( "            S:::::SN::::::N     N:::::::::N    A:::::AAAAAAAAAAAAA:::::A      F:::::F              U::::::U   U::::::U ")
 print( "SSSSSSS     S:::::SN::::::N      N::::::::N   A:::::A             A:::::A   FF:::::::FF            U:::::::UUU:::::::U ")
 print( "S::::::SSSSSS:::::SN::::::N       N:::::::N  A:::::A               A:::::A  F::::::::FF             UU:::::::::::::UU  ")
 print( "S:::::::::::::::SS N::::::N        N::::::N A:::::A                 A:::::A F::::::::FF               UU:::::::::UU    ")
 print( " SSSSSSSSSSSSSSS   NNNNNNNN         NNNNNNNAAAAAAA                   AAAAAAAFFFFFFFFFFF                 UUUUUUUUU      \n")
 
"""
______________$$$$$$$$$$____________________
_____________$$__$_____$$$$$________________
_____________$$_$$__$$____$$$$$$$$__________
____________$$_$$__$$$$$________$$$_________
___________$$_$$__$$__$$_$$$__$$__$$________
___________$$_$$__$__$$__$$$$$$$$__$$_______
____________$$$$$_$$_$$$_$$$$$$$$_$$$_______
_____________$$$$$$$$$$$$$_$$___$_$$$$______
________________$$_$$$______$$$$$_$$$$______
_________________$$$$_______$$$$$___$$$_____
___________________________$$_$$____$$$$____
___________________________$$_$$____$$$$$___
__________________________$$$$$_____$$$$$$__
_________________________$__$$_______$$$$$__
________________________$$$_$$________$$$$$_
________________________$$$___________$$$$$_
_________________$$$$___$$____________$$$$$$
__$$$$$$$$____$$$$$$$$$$_$____________$$$_$$
_$$$$$$$$$$$$$$$______$$$$$$$___$$____$$_$$$
$$________$$$$__________$_$$$___$$$_____$$$$
$$______$$$_____________$$$$$$$$$$$$$$$$$_$$
$$______$$_______________$$_$$$$$$$$$$$$$$$_
$$_____$_$$$$$__________$$$_$$$$$$$$$$$$$$$_
$$___$$$__$$$$$$$$$$$$$$$$$__$$$$$$$$$$$$$__
$$_$$$$_____$$$$$$$$$$$$________$$$$$$__$___
$$$$$$$$$$$$$$_________$$$$$______$$$$$$$___
$$$$_$$$$$______________$$$$$$$$$$$$$$$$____
$$__$$$$_____$$___________$$$$$$$$$$$$$_____
$$_$$$$$$$$$$$$____________$$$$$$$$$$_______
$$_$$$$$$$hg$$$____$$$$$$$$__$$$____________
$$$$__$$$$$$$$$$$$$$$$$$$$$$$$______________
$$_________$$$$$$$$$$$$$$$__________________


"""