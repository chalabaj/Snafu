"""" 
Terachem MPI interface adapted from the ABIN code
Modules: forces_tera.F90 forces_terash.F90
Terachem is very sensitive to the type, lenght and order of transferred data        
We take advantage of NUMPY which can set C-like ordering and data types that are exploited by MPI.Send  
.tobytes for numpy array not necessary with the exception of int8 data type
"""
import traceback
import numpy as np
import os
import sys
import time
from datetime import datetime
sys.path.append('/home/srsen/bin/PYTHON/MPI4PY/mpi4py-3.0.0/build/lib.linux-x86_64-3.6/')
from mpi4py import MPI

try:
    from errors import error_exit
    from default import  (
        max_terachem_time
    )
    from constants import *
except ImportError as ime:
    print("Module {} in inits not found.".format(ime))
    exit(1)

def finish_tera(comm):
    print("Disconnecting Terachem communication.")
    comm.Send(str.encode("0"), dest=0, tag=0 )  # MPI_TAG_EXIT = 0
    # also  MPI_TAG_ERROR = 13 if error occus
    comm.Disconnect()
    return
    
def exit_tera(comm):
    # comm.Disconnect()
    comm.Abort()
    return
 
def tera_connect():
    #call MPI_COMM_CONNECT(port_name, MPI_INFO_NULL, 0, MPI_COMM_SELF, newcomm, ierr)
        
    # this should be moved to init routine
    print("-------TERCHAME INTERFACE------")
    try:
        mpi_port = os.environ['MPI_TERA_PORT']
        print("MPI_TERA_PORT {} found ".format(mpi_port)) 
        comm = MPI.COMM_WORLD.Connect(mpi_port, MPI.INFO_NULL, 0)  #  ,MPI.INFO_NULL     
    except KeyError as PE:
        print("MPI port for communication with TERAPORT was not exported.",
              "Check tera.out.") 
    except Exception as ANYE:
        print("MPI connection to Terachem failed.",
              "\n{}".format(ANYE))
        exit(1)
    else:
        print("Terachem connection established through port: {}.".format(mpi_port))          
    return(comm)

def recieve_tera(comm, natoms, nstates, state, pot_eners, 
                 MO, CiVecs, blob, SMatrix, NAC, TDip, Dip,
                 qmcharges, civec_size, nbf_size, blob_size):
    
    status = MPI.Status()
    
    while not comm.Iprobe(source=MPI.ANY_SOURCE, tag=MPI.ANY_TAG, status=status):
        print("Waiting for Terachem to finish calculations",
              "\n Probe status: {}".format(status.Get_error()))
        time.sleep(1) 
    
    try:
        comm.Recv([pot_eners, nstates, MPI.DOUBLE], source=MPI.ANY_SOURCE,
                  tag=MPI.ANY_TAG, status=status)  
        comm.Recv([TDip, (nstates-1)*3, MPI.DOUBLE], source=MPI.ANY_SOURCE,
                  tag=MPI.ANY_TAG, status=status)  
        comm.Recv([Dip, nstates*3, MPI.DOUBLE], source=MPI.ANY_SOURCE,
                  tag=MPI.ANY_TAG, status=status) 
        comm.Recv([qmcharges, natoms, MPI.DOUBLE], source=MPI.ANY_SOURCE, 
                  tag=MPI.ANY_TAG, status=status)  
        comm.Recv([MO, nbf_size*nbf_size, MPI.DOUBLE], source=MPI.ANY_SOURCE, 
                  tag=MPI.ANY_TAG, status=status)
        comm.Recv([CiVecs, nstates*civec_size, MPI.DOUBLE], source=MPI.ANY_SOURCE,
                  tag=MPI.ANY_TAG, status=status)
        comm.Recv([SMatrix, nstates*nstates, MPI.DOUBLE], source=MPI.ANY_SOURCE,
                  tag=MPI.ANY_TAG, status=status)
        comm.Recv([blob, blob_size, MPI.DOUBLE], source=MPI.ANY_SOURCE,
                  tag=MPI.ANY_TAG, status=status)
        for st1 in range(nstates):
            for st2 in range(st1,nstates):
                comm.Recv([NAC, 3*natoms, MPI.DOUBLE], source=MPI.ANY_SOURCE,
                          tag=MPI.ANY_TAG, status=status)
                if (st1 == st2) and (st1 == state):   
                    xyz = 0
                    for iat in range(0,natoms):
                        fx[iat] = -NAC[xyz]
                        fy[iat] = -NAC[xyz+1]
                        fz[iat] = -NAC[xyz+2]
                        xyz = xyz + 3
                  
    except Exception as excpt:
        print(traceback.format_exc())
        raise RuntimeError("Problem during receiving data from Terachem: {}".format(excpt))
    else:
        print("DATA RECEIVED",
              "\n Probe status: {}".format(status.Get_error()),
              "Energies: {}, Nstates: {}".format(en_array.tolist(), nstates))
        return(fx, fy, fz, pot_eners, MO, CiVecs, blob)
    
def send_tera(comm, natoms, nstates, state, sim_time, x ,y, z
              MO, CiVecs, blob, civec_size, nbf_size, blob_size):
    
    FMSinit = 0
    vels = np.zeros((natoms,3), dtype=np.float64)
    byte_coords = np.dstack((x,y,z))  # x,y,z must stack to single array
    bufints = np.empty(12,order='C',dtype=np.intc)
    bufints[0]=FMSinit
    bufints[1]=natoms
    bufints[2]=1               # doCoup
    bufints[3]=0               # TrajID=0 for SH
    bufints[4]=0               # T_FMS%CentID(1)
    bufints[5]=0               # T_FMS%CentID(2)
    bufints[6]=state  # T_FMS%StateID ! currently not used in fms.cpp
    bufints[7]=0       #  YES if write_wfn write to wfn.bin was don, only for restart 
    bufints[9]=state    # iCalcState-1 ! TC Target State
    bufints[9]=state   # jCalcState-1
    bufints[10]=0              # first_call, not used
    bufints[11]=0              # FMSRestart, not used

    
    #  We need to send only upper-triangle matrix
    #  Diagonal elements are gradients, other NACs
    tocacl = np.zeros((nstates, nstates), dtype=np.intc, order='C')
    uti = np.triu_indices(nstates)   #  upper-triangle indices
    tocacl[state][state] = 1 # gradients for current state only, no NACs

    # Send time
    # Send coordinates
    # Send previous diabatic MOs
    # Send previous CI vecs
    # Only needed for numerical NACME, so send 0 instead 
    # Imaginary velocities for FMS, not needed here, sending zeros...
    try:
        comm.Send([bufints, 12, MPI.INT], 0, 2)
        comm.Send([tocacl[uti], nstates*(nstates-1)/2+nstates, MPI.INT], 0, 2)
        comm.Send([sim_time, 1, MPI.DOUBLE], dest=0, tag=2)
        comm.Send([byte_coords, 3*natoms, MPI.DOUBLE], dest=0, tag=2)
        comm.Send([MO, nbf_size*nbf_size, MPI.DOUBLE], dest=0, tag=2)
        comm.Send([CiVecs, civec_size*nstates, MPI.DOUBLE], dest=0, tag=2)
        comm.Send([blob, blob_size, MPI.DOUBLE], dest=0, tag=2)
        comm.Send([vels, 3*natoms, MPI.DOUBLE], dest=0, tag=2)
        comm.Send([vels , 3*natoms, MPI.DOUBLE], dest=0, tag=2)
        #print(MPI.Status().Get_error())
    except Exception as excpt:
        # any error => RAISE => MPI.ABORT => KILL TERA
        print(traceback.format_exc())
        raise RuntimeError("Problem during sending dat to Terachem: {}".format(excpt))
    else:
        print("Data send.")
        return()

    
def alloc_tera_arrays(civec_size, nbf_size, blob_size, natoms, nstates=4):
    """  ABIN memory allocation of Terachem fields:
          allocate(MO(nbf, nbf)) 
          allocate(MO_old(nbf, nbf))
          allocate(CiVecs(civec,nstate))
          allocate(CiVecs_old(civec,nstate))
          allocate(NAC(natom*3))
          allocate(blob(blobsize))
          allocate(blob_old(blobsize))
          allocate(SMatrix(nstate*nstate))
    """
    MO = np.zeros((nbf_size, nbf_size),dtype=np.float64)
    CiVecs = np.zeros((civec_size, nstates),dtype=np.float64)
    NAC = np.zeros((natoms*3),dtype=np.float64)
    blob = np.zeros((blob_size),dtype=np.float64)
    SMatrix = np.zeros((nstates*nstates),dtype=np.float64)
    qmcharges = np.zeros((natoms),dtype=np.float64)
    TDip =  np.zeros(((nstates-1)*3), dtype=np.float64)
    Dip = np.zeros((nstates*3), dtype=np.float64)
    print("TC arrays allocated.")    
    return(MO, CiVecs, NAC, blob, SMatrix, qmcharges, TDip, Dip)

def move_old2new_terash(MO, MO_old, CiVecs, CiVecs_old, blob, blob_old):
    MO = np.copy(MO_old)
    CiVecs = np.copy(CiVecs_old)
    blob = np.copy(blob_old)
    return(MO, MO_old, CiVecs, CiVecs_old, blob, blob_old)

def tera_init(comm, at_names, natoms, nstates, x,y,z):
    # Initial data transfer to Terachem through MPI  init_terash
    print("Sending initial data to Terachem.")     
    FMSinit = 1
    buffer = np.empty(3,dtype=np.intc) 
    natmm_tera = 0 
    at_names = [(at + " ") if len(at) == 1 else at for at in at_names ]   
    byte_coords = np.dstack((x,y,z))  # x,y,z must stack to single array
    byte_ints = np.array([FMSinit, natoms, natmm_tera],order='C',dtype=np.intc) #.tobytes() #cant be 8 byte     
    byte_natoms = np.array([natoms], dtype=np.int8).tobytes() #must 8 byte number, not int standard -> tobytes           
    byte_names = np.array(at_names, order='C', dtype=np.character) #.tobytes()

    #### ABIN initialize_terache(comm) #####################################
    #  call MPI_Send( natqm, 1, MPI_INTEGER, 0, 2, newcomms(itera), ierr )
    #  call MPI_Send(names, 2*natqm, MPI_CHARACTER, 0, 2, newcomms(itera), ierr )
    #  Send atoms names:  call MPI_Send( names, 2*natqm, MPI_CHARACTER, 0, 2, newcomm, ierr )
    #  Send coordinates: Call MPI_Send( qmcoords, 3*natom, MPI_DOUBLE_PRECISION, 0, 2, newcomm, ierr )
    #  status = MPI.Status()  print("Status: {}, Error: {}".format(status.Get_tag(), status.Get_error()))   
    try:
        print("Sending number of QM atoms.")
        comm.Send( [byte_natoms, 1, MPI.SHORT], dest=0, tag=2 )
        print("Sending QM atom names.") 
        comm.Send( [byte_names, 2*natoms, MPI.CHAR], dest=0, tag=2)
        print("Sending FMS init.")
        comm.Send( [byte_ints, 3, MPI.INT], dest=0, tag=2 ) 
        print("Sending atom names")
        comm.Send([byte_names, 2*natoms, MPI.CHAR], dest=0, tag=2)
        print("Sent initial coordinates to TeraChem.")
        comm.Send([byte_coords, 3*natoms, MPI.DOUBLE], dest=0, tag=2)

        #  Lets wait until Tera finished first ES calc.
        cc = 0
        while not comm.Iprobe(source=MPI.ANY_SOURCE, tag=MPI.ANY_TAG, status=status):
            print("Waiting for Terachem to finish calculations",
                  "\n Probe status {}: {}".format(cc, status.Get_error()))
            time.sleep(1) 
            cc += 1 
            if cc >= max_terachem_time:
                exit_tera(comm)   
                error_exit(15, "Didn't receive data from TC in time during initial comminucation") 
        
        #   civec = np.frombuffer(buffer,dtype=np.intc,count=-1)[0] buffer=bytearray(32*3) 32byte*3fields
        #  .tobytes() or buffer = bytearray(32*3)
        #   Call MPI_Recv( bufints, 3, MPI_INTEGER, MPI_ANY_SOURCE, & MPI_ANY_TAG, newcomm, status, ierr)
        comm.Recv([buffer, 3, MPI.INT],source=MPI.ANY_SOURCE, tag=MPI.ANY_TAG,status=status)
        civec_size = buffer[0]
        nbf_size = buffer[1]
        blob_size = buffer[2]      
        print(civec_size, nbf_size, blob_size)         
    except Exception as excpt:
        # any error => RAISE => MPI.ABORT => KILL TERA
        print(traceback.format_exc())
        exit_tera(comm)
        error_exit(15, str("Error during sending initial TC data {}".format(excpt)))
    else:
        print("TC init done.")
    
    MO, CiVecs, NAC, blob, SMatrix, \
    qmcharges, TDip, Dip = alloc_tera_arrays(civec_size, nbf_size, blob_size,
                                             natoms, nstates)

    return(MO, MO_old, CiVecs, CiVecs_old, NAC, blob, blob_old, SMatrix, 
           civec_size, nbf_size, blob_size, qmcharges, TDip, Dip)
    
    """ 
        if __name__ == "__main__":
        #at_names = ["O","O","H","H ","H ","H "]  # taken from input
        natoms = 6
        nstates = 4
        state = 0
        sim_time = 0.0005
        sim_time = np.array(sim_time, dtype=np.float64)
        # These will be in INIT ROUTINE
        en_array = np.zeros(nstates, dtype=np.float64) 
    
        # np.stack((a, b), axis=-1)
        byte_coords = np.loadtxt("movie.xyz", usecols=(1,2,3),dtype=np.float64)*ANG_BOHR
        print(byte_coords)
    
        comm = tera_connect()
    
        MO, MO_old, CiVecs, CiVecs_old, NAC, blob, \
        blob_old, SMatrix, civec_size, nbf_size,  \
        blob_size, qmcharges, TDip, Dip = tera_init(comm, at_names, natoms, nstates,
                                         byte_coords)
        for i in range(4):
            print("#######{}######".format(i))
            byte_coords = byte_coords *0.99
            sim_time = sim_time * 0.95
            try:
                send_tera(comm, natoms, nstates, state, sim_time, byte_coords, MO,
                          CiVecs, blob, civec_size, nbf_size, blob_size) 
    
                en_array, MO, CiVecs, blob, \
                SMatrix, NAC = recieve_tera(comm, natoms, nstates, en_array, MO,
                                            CiVecs, blob, SMatrix, NAC, TDip, Dip,
                                            qmcharges, civec_size, nbf_size, 
                                            blob_size, state)
                                          
                print(en_array) 
                print(MO)
                print(CiVecs)
                print(SMatrix)
                print(NAC)
                print("\n \n \n ----------------------------------------------")
            except Exception as excpt:
                print("Something went wrong during MPI SEND/RECEIVE.",
                      "\n{}".format(excpt))
                exit_tera(comm)
                exit(1)
        print("ALL DONE OK")
        finish_tera(comm)
        exit(0)
    
    """
  
def calc_forces_tera(comm, natoms, nstates, state, sim_time, x, y, z,
                     fx, fy, fz, pot_eners
                     MO, CiVecs, blob, NAC, civec_size, nbf_size, blob_size)
    try:
        send_tera(comm, natoms, nstates, state, sim_time, x ,y, z
                  MO, CiVecs, blob, civec_size, nbf_size, blob_size)                 
                         
        pot_eners, MO, CiVecs, blob, \
        fx, fy, fz, pot_eners, MO, CiVecs, blob = recieve_tera(comm, natoms, nstates, state, pot_eners, 
                                             MO, CiVecs, blob, SMatrix, NAC, TDip, Dip,
                                             qmcharges, civec_size, nbf_size, blob_size)                 
    except Exception as excpt:
                print("Something went wrong during MPI SEND/RECEIVE.",
                      "\n{}".format(excpt))
                exit_tera(comm)
                error_exit(error_exit(15, str("Error during sending/receive TC data {}".format(excpt))))                             
    return(pot_eners, MO, CiVecs, blob, fx_new, fy_new, fz_new)