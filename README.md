## SNAFU 

### **ab initio** nonadiabatic MD without calculating non-adiabatic couplings

### Theory

1) adiabatic-potential-based formula derived within Landau-Zener:

AK Belyaev, PHYSICAL REVIEW A 84, 014701 (2011), **doi:10.1103/PhysRevA.84.014701**

AK Belyaev, The Journal of Chemical Physics 147, 234301 (2017); **doi:10.1063/1.5000718**


## Status
The code is in **development** stage with no guarantees of its accuracy, precision or fitness. 

* Energy conservation from tests is about 10^-4 - 10^-1 eV between hops  and 10^-5 eV for the regions without hops. 

* Velocity adjustment after hops seems to be more stable when the velocities are scaled by the simple factor K = sqrt(1+-dE/Ekin), where dE is difference between potential energies, for which the hop occured, and Ekin is the kinetic energy at the moment of a hop. Another option is to apply new forces of the final state after hop, however, this requires extra calculations of the forces and the energy conservation appears to be less stable.

* Timestep of 4 au appears to be most suitable for hopping algorithm, but that depend on the PES complexity and some testing is always recommended as to minimize the number of hops. Timestep between 2-8 au should be sufficient.

* Hops are only allowed for pot. energy differences with less than 0.5 eV by default. Higher energy transitions generally lead to poor energy conservation.

## How to run
1) System requirements:
The code was tested on Linux/Debian 4.7.2-5 platform and Anaconda 3.6 package. MPI parallel environment for TeraChem-1.9 point-to-point communication was tested with mpich-3.1.3 and mpi4py-3.0.0. Mpi4py module has to be build with the same mpich as Terachem. 

1) **ABINITIO** folder has to contain one of the .sh script from INTERFACES folder or a terachem input file (or you can add your own interface, see below) .  
Option **abinitio** in input.in must equal to the name of the interface script (e.g. abinitio = molpro-casscf.sh). Interface file script ABINITIO/molpro-casscf.sh must exist. Same apply to terachem input file.
Current implemented interfaces:  
CASSCF in MOLPRO (2015.1)  
CASPT2 in MOLPRO 
TDDFT  in ORCA(4.0.1) & GAUSSIAN 09 (thresh needs to be set otherwise some tddft vectors might not fully converge and still gaussian exits with 0, one needs to check this on the run)
OM3    in MNDO99 
FOMO-CASCI in TERACHEM (works through MPI interface through FMS interface (see TC manual))
BOMD: in Gaussian09 (MP2/DFT) or MOLPRO(MP2)

TODO: EOM-IP, EOM-EA in QCHEM/ORCA,  
** Adding new ab initio interface **
It is straight forward to implement a new ab initio interface as the code reads the gradients.dat file at each step (see examples in INTERFACE folder). This file is created by an interface script after an ab initio code completes energy and gradient calculations. Interface script greps energies and gradients to the gradients.dat file in a running directory with the following structure:  
energy-gs  
energy-1ex state  
....  
energy-nstate  
grad_x(1at) grad_y(1at) grad_z(1at)  
....  
grad_x(n_at) grad_y(n_at) grad_z(n_at)  

Be carefull: the code expects gradients to be extracted, not forces (e.g. Gaussian code). If your ab initio code print forces, you can either grep the negative values of forces to get the gradients (F = - grad(E)) or rename your interface script so that the script name contains one the words forces,gaus or g09 ,separated by dot or minus sign (e.g. forces.sh, forces-gaus.sh, forces-abinitiocode.sh etc), SNAFU will then transform forces into gradients forces.
Also: as the SNAFU is PYTHON based, so the states starts from 0 (ground state). Hence, if an ab initio code denotes ground state as 1, there should be modification in the interface skript (see e.g. molpro examples).  

3) In order to run the code, **geom.in** with the initial geometry and **input.in** files have to present in a running folder.
The **veloc.in*** file with initial velocities is not mandatory. If there are no initial velocities, a simulation will start with zero velocities.
The geom.in and veloc.in files have to be in the XYZ format. Positions are in **ANGSTROM** units (so you can modify it with e.g. molden), but velocities have to be in the **ATOMIC** units:  
Examples of water molecule in geom.in:
3  
water molecule in angstorm units
 O     0.000000     0.000000     0.000000
 H     0.000000     0.000000     0.947000
 H     0.892841     0.000000    -0.315663  

4) Launching:
The LAUNCHER folder contains launchSNAFU and SNAFUS bash scripts which will start the simulation. These are customized for Linux cluster-type computers with queuing systems (e.g. SGE). 
The launchSNAFU is not needed when running SNAFU directly without queuing; SNAFUS script must be adjusted as it takes SGE queuing parameters like JOB_ID env variable. 
The launcher:  
- submit the job to que
- export necessary environment variables (SNAFU_DIR, MPI_TERA)
- create folder on scratch  
- start the simulation on scratch
- copy data back to folder

During initial checks, SNAFU requires to find **SNAFU_DIR** environment variable (in python this is os.environ['SNAFU_DIR']). This variable points to the folder with the **snafu.py** file and **snafu** subfolder containing all the modules. The variable is exported in the SNAFUS launcher script and should be modified depending on where you keep the code. The same applies for Terachem/MPI which are also exported here.

SNAFU simulation with TeraChem interface can be launched by:  
</code>launchSNAFU 1 aq-gpu-gtx980 tera </code>

If the launcher is not used, export **SNAFU_DIR** and MPI_TERA:
  <code>
   export SNAFU_DIR="path/to/snafu/dir"
   export MPI_TERA=0  # MPI_TERA=1 if Terachem interface is used
  </code>
and then run the code:
 <code>
 python snafu.py > snafu.out
 </code> 
This will not work for TeraChem jobs, where terachem runs in background and awaits MPI communication with SNAFU (see TERASNAFUS). 
File snafu.out file is copied during restart for backup. If you used other output filename, it will not becopied during restart backup process.  
 
Environment variables for particular ab-initio code are exported in the interfaces scripts and have to be adjusted to your machine environment (see interface scripts in INTERFACE folder).

## How to restart dynamics

You can restart dynamics from the last completed step or from the chosen step XX depending on how often you wrote restart file during the original simulation.  
The **restart.in** file contains all needed information from the last completed simulation step. Similarly, **restart_freq = XX** option sets the interval for writing a restart file (restart_freq = 100=> restart_100.in, restart_200.in,restart_300.in,...). Before restart, do NOT change the **input.in** file with an exception to extend the simulation time by increasing the **maxsteps** option.   

* To restart simuluation from the last completed step, set **restart = 1** and restart.in file must be in executing folder.
* To restart simuluation from XX step, set **restart = XX** and restart_XX.in file must be in executing folder.

If you restart from some step, existing restart files with the same name will be overwritten (e.g. if you restart from 10th step, all restart files after that step will be overwritten)..

During the restart process, all previous output files (i.e. movie.xyz, energies.dat, restart*.in files, state.dat, snafu.out, velocities.xyz and PES.dat) will be copied to the folder named **PREV_RUN${N}** where N depends a on number of previous restarts (PREV_RUN0 folder contains original simulation data). If you redirect simulation output (python snafu.py > snafu.out) to a file other than the snafu.out, it will not be backed-up. 

The output files are opened in the "append" mode. This will ensure the continuation of output files, however, the original files will rather be backed-up. This is important since the restart procedure truncate all the output files (except of the snafu.out output file which start from empty file) after XX step.
## Input.in options:

[Settings]  
natoms  = 3                # number of atoms in system  
nstates = 3                # number of electronic states  
init_state = 2             # initial electronic state, 0 => ground state, 1 => first ex. state  
timestep = 6               # in atomic unit, 1 au = 0.024 fs   
maxsteps = 600             # total number of steps  
method  = lz               # lz (default)/bomd on selected state(no hops allowed)
abinitio  = molpro-casscf.sh  # name of ab initio interface file in the ABINITIO folder
vel_adj = 0                # 0  - simple scaling K = sqrt(1+-dE/Ekin) default, 1- forces from new surface are included into velocity at hop point    
ener_thresh = 1.0          # threshold for max energy drift (in eV)     
hop_thresh = 0.5           # energy threshold for hopping between the states with energy difference less than this (in eV)    
restart = 0                # N - restart from N-th step, restart_N.in must exist
                           # 1 - restart from the last completed step (i.e. restart.in)
                           # 0 - unset but writes restart information
restart_freq = 100         # writes restart_N.in file each N-th step, here N = 100 (100, 200, etc.) (default = 100)  
write_freq = 100           # how often print output (default 10) 


## TODO:
add diabatization scheme: Le Yu, Phys.Chem.Chem.Phys., 2014, 16, 25883; **doi:10.1039/C4CP03498H**  
