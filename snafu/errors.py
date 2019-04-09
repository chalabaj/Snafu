
import sys
import os

def error_exit(error_number, error_desc=" "):
    """
    Exiting SNAFU simulation.
    MPI4PY implementation of MPI is currently not perfect when dealing with errors and can cause deadlocks. Hence MPI checking.
    """

    err = ("0 - File input.in not found in folder.\n",
           "1 - File geom.in not found in folder.\n",
           "2 - Number of atoms in geoms.in and input.in is not consistent. Please check geom.in and input.in and check XYZ format.\n",
           "3 - OS error: {0}".format("Could not create necessarry files.\n"),
           "4 - Ab initio calculations failed.\n",
           "5 - Ab initio interface (r.X X=gauss, molpro etc.) not found. Check ABINITIO folder.\n",
           "6 - Hopping probability larger than 1, something went wrong.\n",
           "7 - Too large energy drift.\n",
           "8 - File {} exists, but the restart option is turned off (restart = 0). Either remove file or change the restart option.\n".format(error_desc),
           "9 - Input varible(s) is not properly set. See input.in.\n{}".format(error_desc),
           "10 - Restart file {} was not found.".format(error_desc),
           "11 - Wrong input parameter.\nError:{}.".format(error_desc),
           "12 - Input variable {} not set.".format(error_desc),
           "13 - Error in TeraChem interface.Error:{}".format(error_desc),
           "14 - Wrong format of geom.in or veloc.in. Should be XYZ format with geometry in Angstrom and velocities in atomic units.\nError:{}.".format(error_desc),
           "15 - {}".format(error_desc),
           "16 - Error when reading data from restart file.",
           "17 - Error when reading restart file. Option \"{}\" not found.".format(error_desc),
           "18 - User called ABORT {}.".format(error_desc),
           "19 - Import error: {}.".format(error_desc)
          )
    print("-------------------------------------------------------------------")
    print(err[error_number])
    print("Program was terminated due to an error!")
    sys.stdout.flush() 
    
    # prevent MPI deadlock without raiseing runtime error which is caught by excepthook, does not create traceback
    tera_mpi = 0
    tera_mpi = int(os.environ['MPI_TERA'])
    if tera_mpi:
         from tera_propagates import global_except_hook
         sys.excepthook = global_except_hook 
         raise RuntimeError() 
    sys.exit(1)    
    return()
