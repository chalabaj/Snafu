## SNAFU 

### **ab initio** nonadiabatic MD without calculating non-adiabatic couplings

### Theory

1) adiabatic-potential-based formula derived within Landau-Zener:

AK Belyaev, PHYSICAL REVIEW A 84, 014701 (2011), **doi:10.1103/PhysRevA.84.014701**

AK Belyaev, The Journal of Chemical Physics 147, 234301 (2017); **doi:10.1063/1.5000718**


## Status
The code is in **development** stage with no guarantees of its accuracy, precision or fitness. 

* Energy conservation from tests is about 10^-4 eV - 10^-1 eV between hops  and 10^-5 eV for the regions without hops. The WF electronic structure is as the critical point in this context. The CASSCF might produce energy jumps in regions with significant structural changes, so one should examine the PES of a given trajectory. 

* Velocity adjustment after hops seems to be more stable when the velocities are scaled by the simple factor K = sqrt(1+-dE/Ekin), where dE is difference between potential energies, for which the hop occured, and Ekin is the kinetic energy at the moment of a hop. Another option is to apply new forces of the final state after hop, however, this requires extra calculations of the forces and the energy conservation appears to be less stable.

* Timestep of 4 au appears to be most suitable for hopping algorithm, but that depend on the PES complexity and some testing is always recommended as to minimize the number of hops. Timestep between 2-8 au should be sufficient.

* Hops are only allowed for pot. energy differences with less than 0.5 eV by default. Higher energy transitions lead to poor energy conservation.


## How to run
The code was tested on Linux/Debian 4.7.2-5 platform with Anaconda 3.6 package.
Python ver 3.5 and newer should work.

1) Set up environment variable **SNAFU_DIR** in .bashrc or export it before launching the code (see LAUNCHER folder). 
  <code>
   export SNAFU_DIR="path/to/snafu/dir"
  </code>



2) **ABINITIO** folder have to contain one of the script from INTERFACES folder (or you can add your own interface, see below).  
Option **abinitio** in input.in must equal to the name of the script without the .sh extention (e.g. abinitio=molpro-casscf for molpro-casscf.sh script in /ABINITIO folder)
Currently fully working:  
CASSCF in MOLPRO (2015.1)  
ORCA(4.0.1): tddft -  working  
GAUSSIAN 09: tdddft - working (thresh needs to be set otherwise some tddft vectors might not fully converge and still gaussian exits with 0 - weird)
BOMD:  MP2 in gaussian or MOLPRO  
TERACHEM: works through FMS interface, testing


TODO: EOM-IP, EOM-EA in QCHEM/ORCA,  
It is straight forward to implement a new ab initio interface as at each step, the code reads the gradients.dat file in the running directory with the following structure:  
energy-gs  
energy-1ex state  
....  
energy-nstate  
Fx(1at) Fy(1at) Fz(1at)  
....  
Fx(n_at) Fy(n_at) Fz(n_at)  

Environment variables for particular ab initio code have to be adjusted to your machine environment (e.g. if Molpro is used then $MOLPROEXE variable should be set up).


3) In order to run the code, **geom.in** with the initial geometry and **input.in** files have to present in a folder.
The **veloc.in*** file with initial velocities can be also used, otherwise the dynamics will start with zero velocities.
The geom.in and veloc.in files have to be in the XYZ format and velocities in atomic units:  
3  
water xyz  
O 0.0 0.0 0.0   
H 1.0 0.0 0.0   
H 0.0 1.0 0.0   


## How to restart dynamics

You can restart dynamics from the last completed step or from the chosen step depending on how often you wrote restart file in an original simulation.  
The restart.in file contains all needed information from the last completed simulation step. Similarly, checkpoints are created according to the **restart_freq** option which sets the interval for writing a restart file (restart_400.in contains restart information from the 400th step). File **input.in** should not be changed with an exception to extend the simulation time by increasing the **maxsteps** option.

* To restart simuluation from the last completed step, set **restart = 1** and restart.in file must be in executing folder.
* To restart simuluation from XX step, set **restart = XX** and restart_XX.in file must be in executing folder.

During each restart, all previous output files (i.e. movie.xyz, energies.dat, restart*.in, input.in, state.dat, snafu.out, velocities.xyz and PES.dat) will be copied to the folder named **PREV_RUN${N}** where N depends on number of previous restarts (PREV_RUN0 folder contains original simulation data).
Copying files is executed on the launcher (launchSNAFU) level, is if one uses other launcher or runs the code directly with queing system, there will be no back-up of original files and you can LOST trimmed data (if restart > 1).

The output files are opened in the "append" mode. This will ensure the continuation of output files, however, the original files will rather be backed-up. This is important since the restart procedure trims the output files (except of the snafu.out) after XX step.
## Input.in options:

[Settings]  
natoms  = 3                # number of atoms in system  
nstates = 3                # number of electronic states  
init_state = 2             # initial electronic state, 0 => ground state, 1 => first ex. state  
timestep = 6               # in atomic unit au = 0.024 fs   
maxsteps = 600             # total number of steps  
method  = lz-adiabatic     # bomd/lz-adibatic (Belyaev)
abinitio  = molpro-casscf.sh  # ab initio interface file, has to start:  g09, molpro, orca, tera input file (e.g. tera.inp)
vel_adj = 1                # 0  - simple scaling K = sqrt(1+-dE/Ekin), 1- forces from new surface are included into velocity at hop point    
ener_thresh = 1.0          # threshold for max energy drift in eV     
hop_thresh = 0.5           # energy threshold for hopping between the states with energy difference less than this (in eV)    
restart = 0                # N - restart from N-th step, restart_N.in must exist
                           # 1 - restart from the last completed step (i.e. restart.in)
                           # 0 - unset but writes restart information
restart_freq = 100         # writes restart_N.in file each N-th step, here N = 100 (100, 200, 300 etc.),default = 100
tera_mpi = 0               # use Terachem abinitio interface via MPI   
write_freq = 100           # how often print output, default 10
## TODO:
add diabatization scheme: Le Yu, Phys.Chem.Chem.Phys., 2014, 16, 25883; **doi:10.1039/C4CP03498H**  
add restart option
