# Terachem interface adapted from ABIN code

import numpy as np
import os
import sys
import time
from datetime import datetime
#sys.path.append('/home/srsen/bin/PYTHON/MPI4PY/mpi4py-3.0.0/build/lib.linux-x86_64-3.6/')
start = 1 # first geometry to process
from mpi4py import MPI
end = 0 # last geometry to process. set to 0 to process all geometries to the end of the movie.
name = "NH3"
suffix = "com.xyz"
movie = "movie.xyz" # file with geometries to process
output = "energies_ground.dat" # name of the output file with resulting energies

def finish_tera(comm):
    print("Exiting MPI Terachem communication.")
    comm.Send(str.encode("0"), dest=0, tag=0 )  # MPI_TAG_EXIT = 0
    # also  MPI_TAG_ERROR = 13 if error occus
    comm.Disconnect()
    return()
 
def tera_connect():
    #call MPI_COMM_CONNECT(port_name, MPI_INFO_NULL, 0, MPI_COMM_SELF, newcomm, ierr)
        
    # this should be moved to init routine
    print("Trying to connect to Terachem server:")
    try:
        mpi_port = os.environ['MPI_TERA_PORT']
        print("MPI_TERA_PORT {} found ".format(mpi_port))                
    except KeyError as PE:
        print("MPI port for communication with TERAPORT was not exported. Check tera.out.")
        
    comm = MPI.COMM_WORLD.Connect(mpi_port, MPI.INFO_NULL, 0)  #  ,MPI.INFO_NULL
    if comm:
        print("Terachem connection established.")
    else:
        print("MPI connection to Terachem failed.")
        exit(1)
    return(comm)

def alloc_tera_arrays(civec, nbf, blobsize, natoms, nstates=4):
    
    startTime = datetime.now()

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
    MO = np.empty((nbf, nbf),dtype=np.float64)
    MO_old = np.empty((nbf, nbf),dtype=np.float64)
    CiVecs = np.empty((civec,nstates),dtype=np.float64)
    CiVecs_old = np.empty((civec,nstates),dtype=np.float64)
    NAC = np.empty((natoms*3),dtype=np.float64)
    blob = np.empty((blobsize),dtype=np.float64)
    blob_old = np.empty((blobsize),dtype=np.float64)
    SMatrix = np.empty((nstates*nstates),dtype=np.float64)
    stopTime = datetime.now()
    simtime = (datetime.now() - startTime)
    print("Tera arrays allocation times: (hh:mm:ss): {}".format(simtime))

    return(MO, MO_old, CiVecs, CiVecs_old, NAC, blob, blob_old, SMatrix)

def move_old2new_terash(MO, MO_old, CIVecs, CIVecs_old, blob, blob_old):
    MO = np.copy(MO_old)
    CIVecs = np.copy(CIVecs_old)
    blob = np.copy(blob_old)
    return(MO, MO_old, CIVecs, CIVecs_old, blob, blob_old)

def tera_init(comm, natoms = 6):
    """ Initial data transfer to Terachem through MPI
    Terachem is very sensitive to type, lenght and order of transferred data        
    We take advantage of NUMPY which can set C-like ordering and data types
    .tobytes for numpy arraay is not needed as numpy keep daty in byte anyways, yet dtype must be set     
    """   
    print("Sending data to Terachem:")        
    FMSinit = 1
    natoms = 6
    natmm_tera = 0 
    at_names = ["O","O","H","H ","H ","H "]  # taken from input
    at_names = [(at + " ") if len(at) == 1 else at for at in at_names ]   
   
    byte_ints = np.array([FMSinit, natoms, natmm_tera],order='C',dtype=np.intc).tobytes() #cant be 8 byte     
    byte_natoms = np.array([natoms], dtype=np.int8).tobytes() #must be 8 byte number...           
    byte_names = np.array(at_names, order='C', dtype=np.character) #.tobytes()

    #### ABIN initialize_terache(comm) #####################################
    print("Sending initial number of QM atoms to TeraChem.")
    #  call MPI_Send( natqm, 1, MPI_INTEGER, 0, 2, newcomms(itera), ierr )
    comm.Send( [byte_natoms, 1, MPI.SHORT], dest=0, tag=2 )
    
    print("Sending initial QM atom names to TeraChem.") 
    #  call MPI_Send(names, 2*natqm, MPI_CHARACTER, 0, 2, newcomms(itera), ierr )
    comm.Send( [byte_names, 2*natoms, MPI.CHAR], dest=0, tag=2)

    # ABIN  init_terash ####################################################
    comm.Send( [byte_ints, 3, MPI.INT], dest=0, tag=2 ) 
    print("Sent initial FMSinit.")
    status = MPI.Status()
    if status.Get_error():
        finish_tera(comm)
        exit(1)
    print("Status: {}, Error: {}".format(status.Get_tag(), status.Get_error()))
   
    #  Send atoms names
    #  call MPI_Send( names, 2*natqm, MPI_CHARACTER, 0, 2, newcomm, ierr )
    comm.Send([byte_names, 2*natoms, MPI.CHAR], dest=0, tag=2)
    print("Sending QM atom names to TeraChem.")

    #  Send coordinates
    byte_coords = np.loadtxt(movie, usecols=(1,2,3),dtype=np.float64)  #.tobytes('C')
    # call MPI_Send( qmcoords, 3*natom, MPI_DOUBLE_PRECISION, 0, 2, newcomm, ierr )
    comm.Send([byte_coords, 3*natoms, MPI.DOUBLE], dest=0, tag=2)
    print("Sent initial coordinates to TeraChem.")
    status = MPI.Status()
    print("Status: {}, Error: {}".format(status.Get_tag(), status.Get_error()))
    
    #call MPI_Recv( bufints, 3, MPI_INTEGER, MPI_ANY_SOURCE, &
    #               MPI_ANY_TAG, newcomm, status, ierr)
    while not comm.Iprobe(source=MPI.ANY_SOURCE, tag=MPI.ANY_TAG):
        print("Waiting for Terachem to finish calculations")
        time.sleep(1)            
    buffer = bytearray(32*3)
    buffer = np.empty(3,dtype=np.intc) #.tobytes()
    comm.Recv([buffer, 3, MPI.INT],source=MPI.ANY_SOURCE, tag=MPI.ANY_TAG,status=status)
    #civec = np.frombuffer(buffer,dtype=np.intc,count=-1)[0]
    #nbf   = np.frombuffer(buffer,dtype=np.intc,count=-1)[1]
    #blobsize = np.frombuffer(buffer,dtype=np.intc,count=-1)[2]
    civec = buffer[0]
    nbf = buffer[1]
    blobsize = buffer[2]      
    print(civec,nbf, blobsize)
    
    MO, MO_old, CiVecs, CiVecs_old, NAC, blob, blob_old, SMatrix = alloc_tera_arrays(civec, nbf, blobsize, natoms, nstates = 4)
    blob = 0.0
    blob_old = 0.0
    return()
    
if __name__ == "__main__":
    comm = tera_connect()
    tera_init(comm)
    finish_tera(comm)
    exit(0)