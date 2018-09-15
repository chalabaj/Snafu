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
Energy conservation from initial tests is about 10^-2 eV between hops and 10^5 eV for bomd without hops.
Velocity adjustment after hops seems to be more stable when the velocities are scales factor K = sqrt(1+-dE/Ekin), where dE is difference between potential energies for which the hop occured and Ekin is current kinetic energy at the moment of the hop. Another option is to apply new forces of the final state after hop, however, this requires extra calculations of the forces and the energy conservation is  only about 10^-1 eV.
Timestep of 4 au appers to be most stable, but that depend on the PES complexity and some testing is always recommended as to minimize the number of hops.
Hops are only allowed for pot. eneergy differences with less than 0.5 eV. Higher energy transitions lead to poore energy conservation.


## How to run/Requirements
--very provisionary 
The code was tested on Linux/Debian 4.7.2-5 platform with Anaconda 3.6 package.
Python ver 3.5 and newer should work.

1) Set up environmentn your environment
```
export SNAFU_DIR="path/to/snafu/dir"
``` 

2) ABINITIO folder should contain one of the script from INTERFACES folder, depending on the ab initio code you use:
Currently fully working is CASSCF in MOLPRO
BOMD with MP2 in Gaussian
ORCA (tddft - not fully working yet.
Environment variables for particular ab initio code have to be adjusted to your machine environment (e.g. if Molpro is used then $MOLPROEXE variable should be set up).

## Input variables
File names input.in needs to be present in the folder with following parameters:

[Settings]
natoms  = 3                # number of atoms in system
nstates = 3                # number of electronic states
init_state = 2             # initial electronic state, 0 => ground state, 1 => first ex. state
timestep = 6               # in atomic unit au = 0.024 fs 
maxsteps = 600             # total number of steps
method  = lz-adiabatic     # bomd/lz-adibatic (Belyaev)
abinitio  = molpro-casscf   # where to take gradients, file has to start  g09 or molpro to distinguish between forces and gradients
restart = 0                # not yet working
vel_adj = 1                # 0  - simple scaling K = sqrt(1+-dE/Ekin), 1- forces fro new surface are included into velocity at hop point
ener_thresh = 1.0          # threshold for max energy drift in eV 
