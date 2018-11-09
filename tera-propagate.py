# Terachem interface adapted from ABIN code

import numpy as np
import os
import sys

#sys.path.append('/home/srsen/bin/PYTHON/MPI4PY/mpi4py-3.0.0/build/lib.linux-x86_64-3.6/')
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
    #call MPI_COMM_CONNECT(port_name, MPI_INFO_NULL, 0, MPI_COMM_SELF, newcomm, ierr)
        
    # this should be moved to init routine
    try:
        mpi_port = os.environ['MPI_TERA_PORT']
    except KeyError as PE:
        print("MPI port for communication with TERAPORT was not exported. Check tera.out.")
    
    comm = MPI.COMM_WORLD.Connect(mpi_port, MPI.INFO_NULL, 0)  #  ,MPI.INFO_NULL
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
    bufint = np.array([FMSinit, natoms, natmm_tera],dtype=np.intc) # cant be 8 byte
    byte_atom = np.array([natoms],dtype=np.int8).tobytes()  # must be 8 byte number...
    byte_ints = bufint.tobytes()
    
    
    atoms = ""
    with open(movie, 'r') as file:
	      for i in range(0,natoms):
		        atoms += file.readline()[0:2]    
    file.close()
    byte_names = str.encode(atoms)
    
    #atom_names = np.array(["O","O ","H ","H ","H ","H "],order='C',)
    #byte_names = atom_names.tobytes()
    
    # subroutine  initialize_terache(comm):
    print("Sending initial number of QM atoms to TeraChem.")
    #  call MPI_Send( natqm, 1, MPI_INTEGER, 0, 2, newcomms(itera), ierr )
    comm.Send( byte_atom, dest=0, tag=2 )
    
    print("Sending initial QM atom names to TeraChem.") 
    #  call MPI_Send(names, 2*natqm, MPI_CHARACTER, 0, 2, newcomms(itera), ierr )
    comm.Send( [byte_names, 2*natoms, MPI.CHARACTER], dest=0, tag=2)

    # subroutine  init_terash:
    
    comm.Send( [byte_ints, 3, MPI.INTEGER], dest=0, tag=2 ) 
    print("Sent initial FMSinit. Status: {}".format(MPI.Status()))
    status = MPI.Status()
    if status.Get_error():
        finish_tera(comm)
        exit(1)
    print("Status: {}, Error: {}".format(status.Get_tag(), status.Get_error()))
   
   
    #  Send atoms names
    #  call MPI_Send( names, 2*natqm, MPI_CHARACTER, 0, 2, newcomm, ierr )
    comm.Send([byte_names, 2*natoms, MPI.CHARACTER], dest=0, tag=2)
    print("Sending initial QM atom names to TeraChem.")
    status = MPI.Status()
    print("Status: {}, Error: {}".format(status.Get_tag(), status.Get_error()))

    print("Status: {}, Error: {}".format(status.Get_tag(), status.Get_error()))
    #  Send coordinates
    coords = np.loadtxt(movie, usecols=(1,2,3),dtype=np.float64)
    print(coords)
    # call MPI_Send( qmcoords, 3*natom, MPI_DOUBLE_PRECISION, 0, 2, newcomm, ierr )
    byte_coords = coords[0:natoms, :].tobytes('C')
    print(byte_coords)
    comm.Send([coords, 3*natoms, MPI.DOUBLE], dest=0, tag=2)
    print("Sent initial coordinates to TeraChem.")
    status = MPI.Status()
    print("Status: {}, Error: {}".format(status.Get_tag(), status.Get_error()))
    

    #call MPI_Recv( bufints, 3, MPI_INTEGER, MPI_ANY_SOURCE, &
    #               MPI_ANY_TAG, newcomm, status, ierr)
    buffer = bytearray(8*4*3)
    comm.Recv([buffer, 3, MPI.INTEGER],source=MPI.ANY_SOURCE, tag=MPI.ANY_TAG,status=status)
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