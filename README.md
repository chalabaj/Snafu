## SNAFU 

### **ab initio** nonadiabatic MD without calculating non-adiabatic couplings

### Theory

1) adiabatic-potential-based formula derived within Landau-Zener:

AK Belyaev, PHYSICAL REVIEW A 84, 014701 (2011), **doi:10.1103/PhysRevA.84.014701**

AK Belyaev, The Journal of Chemical Physics 147, 234301 (2017); **doi:10.1063/1.5000718**


## Status
The code is in **development** stage with no guarantees of its accuracy, precision or fitness. 

* Energy conservation from tests is about 10^-4 eV - 10^-1 eV between hops  and 10^5 eV for the regions without hops. The WF electronic structure is as the critical point in this context. The CASSCF might produce energy jumps in regions with significant structural changes, so one should examine the PES of a given trajectory. 

* Velocity adjustment after hops seems to be more stable when the velocities are scaled by the simple factor K = sqrt(1+-dE/Ekin), where dE is difference between potential energies, for which the hop occured, and Ekin is the kinetic energy at the moment of a hop. Another option is to apply new forces of the final state after hop, however, this requires extra calculations of the forces and the energy conservation appears to be less stable.

* Timestep of 4 au appears to be most suitable for hopping algorithm, but that depend on the PES complexity and some testing is always recommended as to minimize the number of hops. Timestep between 2-8 au should be sufficient.

* Hops are only allowed for pot. energy differences with less than 0.5 eV by default. Higher energy transitions lead to poor energy conservation.


## How to run
The code was tested on Linux/Debian 4.7.2-5 platform with Anaconda 3.6 package.
Python ver 3.5 and newer should work.

1) Set up environment, i.e.  set an envir. variable SNAFU_DIR in .bashrc or before launching the code by:
<code>
export SNAFU_DIR="path/to/snafu/dir"
</code>

## How to restart dynamics

* Restart option must be set to 1 and restart.in file must be present in an executing folder.
The restart.in file contains all needed information from the last completed simulation step.   
* If one wants to restart dynamics from other step in history (not the last one), one must use option restart_step (e.g. restart_step = 500) and program will search for file restart_500.in 

2) ABINITIO folder should contain one of the script from INTERFACES folder, depending on the ab initio code you use:  
Currently fully working:  
CASSCF in MOLPRO (2015.1)  
ORCA(4.0.1): tddft -  working  
GAUSSIAN 09: tdddft - working (thresh needs to be set otherwise some tddft vectors might not fully converge  
BOMD:  MP2 in gaussian or MOLPRO  

TODO: EOM-IP, EOM-EA in QCHEM/ORCA, 
It is straight forward to implement a new ab initio as at each step, the code reads gradients.dat file in the running directory with a structure:  
energy-gs  
energy-1ex state  
....  
energy-nstate  
Fx(1at) Fy(1at) Fz(1at)  
....  
Fx(n_at) Fy(n_at) Fz(n_at)  


Environment variables for particular ab initio code have to be adjusted to your machine environment (e.g. if Molpro is used then $MOLPROEXE variable should be set up).

3) In order to run the code, **geom.in** with the initial geometry and **input.in** files have to present in a folder.
The veloc.in files with initial velocities can be also used, otherwise the dynamics will start with zero velocities.
The geom.in and veloc.in files have to be in the XYZ format:  
3  
water xyz  
O 0.0 0.0 0.0   
H 1.0 0.0 0.0   
H 0.0 1.0 0.0   


## Input.in options:

[Settings]  
natoms  = 3                # number of atoms in system  
nstates = 3                # number of electronic states  
init_state = 2             # initial electronic state, 0 => ground state, 1 => first ex. state  
timestep = 6               # in atomic unit au = 0.024 fs   
maxsteps = 600             # total number of steps  
method  = lz-adiabatic     # bomd/lz-adibatic (Belyaev)
abinitio  = molpro-casscf   # where to take gradients and energies, file name has to start  g09 or molpro
vel_adj = 1                # 0  - simple scaling K = sqrt(1+-dE/Ekin), 1- forces from new surface are included into velocity at hop point    
ener_thresh = 1.0          # threshold for max energy drift in eV     
hop_thresh = 0.5           # energy threshold for hopping between the states with energy difference less than this (in eV)    
restart = 0                # 1 - restart from restart.in, 0 - unset but writes restart information
restart_write = 100        # writes restart_$i.in file each i-th step, here i = 100
restart_step = 500         # restart a simulation from selected step, (here 500-th step, restart_500.in must be in executing folder)

## TODO:
add diabatization scheme: Le Yu, Phys.Chem.Chem.Phys., 2014, 16, 25883; **doi:10.1039/C4CP03498H**  
add restart option
