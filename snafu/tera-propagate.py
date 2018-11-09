# Python 3 program for calculating bunch of jobs fast using terachem server
# Stepan Srsen

import numpy as np
import os
import sys

sys.path.append('/home/srsen/bin/PYTHON/MPI4PY/mpi4py-3.0.0/build/lib.linux-x86_64-3.6/')
start = 1 # first geometry to process
from mpi4py import MPI
end = 0 # last geometry to process. set to 0 to process all geometries to the end of the movie.
name = "NH3"
suffix = "com.xyz"
movie = "movie.xyz" # file with geometries to process
output = "energies_ground.dat" # name of the output file with resulting energies

def finish_tera(comm):
    comm.Disconnect()
    return()
 
def tera_connect():
    # forces_tera onnect:
    #call MPI_COMM_CONNECT(port_name, MPI_INFO_NULL, 0, MPI_COMM_SELF, newcomm, ierr)
        
    # this should be moved to init routine
    try:
        mpi_port = os.environ['MPI_TERA_PORT']
    except KeyError as PE:
        print("MPI port for communication with TERAPORT was not exported. Check tera.out.")
    
    comm = MPI.COMM_SELF.Connect(mpi_port, MPI.INFO_NULL, 0)  #  ,MPI.INFO_NULL
    if comm:
        print("Terachem connection established.")
    else:
        print("MPI connection to Terachem failed.")
        exit(1)
    return(comm)

def tera_init(comm):
    FMSinit = 1
    natoms = 6
    natmm_tera = 0
    bufint = np.array([FMSinit, natoms-natmm_tera, natmm_tera],dtype=np.int32)
    bufints = bufint.tobytes()
    bufints=bytearray([FMSinit, natoms-natmm_tera, natmm_tera])
    comm.Send(bufints, dest=0,tag=2) 
    print("Sent initial FMSinit. Status: {}".format(MPI.Status()))
    status = MPI.Status()
    print("Status: {}, Error: {}".format(status.Get_tag(), status.Get_error()))
    
    #  Send atoms names    
    atoms = np.array(["O", "O", "H", "H", "H", "H"],dtype=np.character)
    batoms = atoms.tobytes()
    comm.Send(batoms, dest=0,tag=2)
    print("Sending initial QM atom names to TeraChem.")
    status = MPI.Status()
    print("Status: {}, Error: {}".format(status.Get_tag(), status.Get_error()))
    
    #  Send coordinates
    coords = np.loadtxt(movie, usecols=(1,2,3))
    # call MPI_Send( qmcoords, 3*natom, MPI_DOUBLE_PRECISION, 0, 2, newcomm, ierr )
    comm.Send(coords.tobytes(),dest=0,tag=2)
    print("Sent initial coordinates to TeraChem.")
    status = MPI.Status()
    print("Status: {}, Error: {}".format(status.Get_tag(), status.Get_error()))
    

    #call MPI_Recv( bufints, 3, MPI_INTEGER, MPI_ANY_SOURCE, &
    #               MPI_ANY_TAG, newcomm, status, ierr)
    buffer = bytearray(8*3*2)
    comm.Recv([buffer,MPI.INT],source=0,status=status)
    civec = buffer[1]
    nbf = buffer[2]
    blobsize = buffer[3]
    print(civec,nbf, blobsize)
   
    finish_tera(comm)
    return()
    
if __name__ == "__main__":
    comm = tera_connect()
    tera_init(comm)
    finish_tera(comm)
    exit(0)