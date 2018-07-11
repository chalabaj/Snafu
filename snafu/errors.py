
import sys
import os

def error_exit(error_number):
    err = ("0 - File input.in not found in folder.",
           "1 - File geom.in not found in folder.\n",
           "2 - Number of atoms in geoms.in and input.in is not consistent. Please check geom.in and input.in and check XYZ format.\n",
           "3 - OS error: {0}".format("Could not create necessarry files. Exiting..."))
    print("---------------------------------")
    print("Program was terminated due to an error!")
    print(err[error_number])
    sys.exit(1)
    return()