
import sys
import os

def error_exit(error_number):
    
    """
    TO DO organize error according to time from start
    0 should be init error
    1 should be propagation error...
    """
    err = ("0 - File input.in not found in folder.",
           "1 - File geom.in not found in folder.\n",
           "2 - Number of atoms in geoms.in and input.in is not consistent. Please check geom.in and input.in and check XYZ format.\n",
           "3 - OS error: {0}".format("Could not create necessarry files."),
           "4 - Ab initio calculations failed.",
           "5 - Ab initio interface (r.X X=gauss, molpro etc.) not found. Check ABINITIO folder.")
    print("---------------------------------")
    print("Program was terminated due to an error! Exiting")
    print(err[error_number])
    sys.exit(1)
    return()