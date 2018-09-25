
import sys
import os

def error_exit(error_number):
    
    """
    TO DO organize error according to time from start
    0 should be init error
    1 should be propagation error...
    """
    err = ("0 - File input.in not found in folder.\n",
           "1 - File geom.in not found in folder.\n",
           "2 - Number of atoms in geoms.in and input.in is not consistent. Please check geom.in and input.in and check XYZ format.\n",
           "3 - OS error: {0}".format("Could not create necessarry files.\n"),
           "4 - Ab initio calculations failed.\n",
           "5 - Ab initio interface (r.X X=gauss, molpro etc.) not found. Check ABINITIO folder.\n",
           "6 - Hopping probability larger than 1, something went wrong.\n",
           "7 - Too large energy drift.\n",
           "8 - Some of the output files already exists and restart = 0)\n")
    print("---------------------------------")
    print("Program was terminated due to an error! Exiting...")
    print(err[error_number])
    sys.exit(1)
    return()