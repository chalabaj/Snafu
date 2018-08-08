## SNAFU 

### * *ab initio* * nonadiabatic MD without calculating non-adiabatic couplings

### Theory
Two scheme with/without diabatization

1) adiabatic-potential-based formula derived within Landau-Zener:

AK Belyaev, PHYSICAL REVIEW A 84, 014701 (2011), **doi:10.1103/PhysRevA.84.014701**

AK Belyaev, The Journal of Chemical Physics 147, 234301 (2017); **doi:10.1063/1.5000718**

TODO:
2) diabatization scheme: Le Yu, Phys.Chem.Chem.Phys., 2014, 16, 25883; **doi:10.1039/C4CP03498H**

## Status
Program works in beta version. 

Basic features works, but some needs more testing.
BOMD works with Velocity Verlet.
Hopping algorithm works.
Energy conservation from initial tests is about 10^-2 eV between hops and 10^5 eV between hops, but this need more testing for different ab-initio methods.



## How to run/Requirements
--very provisionary 
The code was tested on Linux/Debian 4.7.2-5 platform with Anaconda 3.6 package.
Python ver 3.5 and newer should work.


1) Set up environmentn your environment
```
export SNAFU_DIR="path/to/snafu/dir"
``` 

2) ABINITIO folder should contain one of the script from INTERFACES folder, depending on the ab initio code you can use (Molpro, Gaussian).
Environment variables for particular ab initio code have to be adjusted to your machine environment (e.g. if Molpro is used then $MOLPROEXE variable should be set up).
