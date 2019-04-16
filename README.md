# SNAFU 

#### Snafu is ab initio surface hopping molecular dynamics without non-adiabatic couplings based on the adiabatic-potential-based formula derived within Landau-Zener approach.

**Theory:** 

AK Belyaev, PHYSICAL REVIEW A 84, 014701 (2011), **doi:10.1103/PhysRevA.84.014701**  
AK Belyaev, The Journal of Chemical Physics 147, 234301 (2017); **doi:10.1063/1.5000718**  


Snafu is python-based code. It was developed and tested on Linux/Debian 4.7.2-5 platform and Anaconda 3.6 package (python 3.6.8). MPI parallel environment for TeraChem-1.9 point-to-point communication was tested with mpich-3.1.3 and mpi4py-3.0.0 package. MPI4PY has to be build with the same mpich version as TeraChem. 

---

## How to run:  

Create a folder containing input files and ABINITIO folder:    

**ABINITIO** folder has to contain one of the .sh script from the INTERFACES folder or a terachem input file (or you can add your own interface, see bellow). Option **abinitio** in input.in must equal to the name of the interface script (e.g. abinitio = molpro-casscf.sh and script ABINITIO/molpro-casscf.sh must exist, same apply to terachem input file).

**Current implemented interfaces:**    
* CASSCF in MOLPRO (2015.1)   
* CASPT2 in MOLPRO   
* TDDFT  in ORCA(4.0.1) & GAUSSIAN 09 (thresh needs to be set otherwise some tddft vectors might not fully converge and still gaussian exits with 0, one needs to check this on the run)  
* OM3    in MNDO99   
* FOMO-CASCI in TERACHEM (works through MPI interface through FMS interface (see TC manual))  
* BOMD: in Gaussian09 (MP2/DFT) or MOLPRO(MP2)  
 
**Input files:**  
In order to run the code, **geom.in** with the initial geometry and **input.in** files have to present in a running folder.
The **veloc.in*** file with initial velocities is not mandatory. If there are no initial velocities, a simulation will start with zero velocities.
The geom.in and veloc.in files have to be in the XYZ format. Positions are in **ANGSTROM** units (so you can modify it with e.g. molden), but velocities have to be in the **ATOMIC** units:  
Examples of water molecule in geom.in:  
```bash
3  
water molecule in angstorm units
 O     0.000000     0.000000     0.000000
 H     0.000000     0.000000     0.947000
 H     0.892841     0.000000    -0.315663  
```

### Input.in options:

```bash
[Settings]  
natoms  = 3                # number of atoms in system  
nstates = 3                # number of electronic states  
init_state = 2             # initial electronic state, 0 => ground state, 1 => first ex. state  
timestep = 6               # in atomic unit, 1 au = 0.024 fs   
maxsteps = 600             # total number of steps  
method  = lz               # lz (default)/bomd on selected state(no hops allowed)
abinitio  = molpro-casscf.sh  # name of ab initio interface file in the ABINITIO folder
vel_adj = 0                # 0  - simple scaling K = sqrt(1+-dE/Ekin) default, 1- forces from new surface are included into velocity after hopping event    
ener_thresh = 1.0          # threshold for max energy drift from start (in eV)     
hop_thresh = 0.5           # energy threshold for hopping between the states with energy difference less than this (in eV)    
restart = 0                # N - restart from N-th step, restart_N.in must exist
                           # 0 - regular run, writes restart information
restart_freq = 100         # writes restart_N.in file each N-th step, here N = 100 (100, 200, etc.) (default = 100)  
write_freq = 100           # how often print output (default 10) 
```
* SNAFU is PYTHON based, so the states starts from 0 (ground state)

* Energy conservation from tests is about 10^-4 - 10^-1 eV between hops  and 10^-5 eV for the regions without hops. 

* Velocity adjustment after hops seems to be more stable when the velocities are scaled by the simple factor K = sqrt(1+-dE/Ekin), where dE is difference between potential energies, for which the hop occured, and Ekin is the kinetic energy at the moment of a hop. Another option is to apply new forces of the final state after hop, however, this requires extra calculations of the forces and the energy conservation appears to be less stable.

* Timestep of 4 au appears to be most suitable for hopping algorithm, but that depend on the PES complexity and some testing is always recommended as to minimize the number of hops. Timestep between 2-8 au should be fine.


### Running  
You can use launchers in the LAUNCHER folder and launchSNAFU and SNAFUS bash scripts which will start the simulation. These are customized for Linux cluster-type computers with queuing systems (here SGE) and should be adjusted to your environment and folder. 
The launcher:  
- submit the job to a que
- export necessary environment variables (SNAFU_DIR, MPI_TERA)
- create folder on scratch  
- start the simulation on scratch
- copy data back to folder


During initial checks, SNAFU requires to find **SNAFU_DIR** environment variable. This variable points to the folder with the **snafu.py** file and **snafu** subfolder containing all the modules. The variable is exported in the SNAFUS launcher script and should be modified depending on where you keep the code. The same applies for Terachem/MPI variable which are also exported here.   
Environment variables for particular ab-initio code are exported in the interfaces scripts and have to be adjusted to your machine environment (see interface scripts in INTERFACE folder).

The launchSNAFU is not needed when running SNAFU directly without que (nodes). If that is the case, export **SNAFU_DIR** and MPI_TERA env variables:  
```bash
   export SNAFU_DIR="path/to/snafu/dir"\
   export MPI_TERA=0  # MPI_TERA=1 if Terachem interface is used  
```  

and then run the code:  
```bash
 python snafu.py > snafu.out
```

This will not work for **TeraChem** jobs, in which the TeraChem runs in background and awaits MPI communication with snafu (see TERASNAFUS). SNAFU simulation with TeraChem interface can be launched by:  
  
```bash
launchSNAFU 1 quename tera
```

---

## How to restart dynamics

You can restart dynamics from the XXth step depending on how often you wrote restart files during the original simulation.  
The **restart_XX.in** file contains all needed information from the XX simulation step. Option **restart_freq = XX** sets the interval for writing a restart file (restart_freq = 100=> restart_100.in, restart_200.in,restart_300.in etc). 

* To restart simuluation from XX step, set **restart = XX** and restart_XX.in file must be in executing folder.

If you restart from some step, existing restart files with the same name will be overwritten (e.g. if you restart from 10th step, all restart files after that step will be overwritten). However, during the restart process, all previous output files (i.e. movie.xyz, energies.dat, restart*.in files, state.dat, snafu.out, velocities.xyz and PES.dat) will be copied to the folder named **PREV_RUN${N}** where N depends a on number of previous restarts (PREV_RUN0 folder contains original simulation data).  

The output files are opened in the "append" mode. This will ensure the continuation of the output files. The restart procedure truncate all the output files after XX step, but the original data are still preserved in a PREV_RUN folder.

---

### Adding new ab initio interface  
It is straight forward to implement a new ab initio interface as the code reads the gradients.dat file at each step (see examples in INTERFACE folder). This file is created by an interface script after an ab initio code completes energy and gradient calculations. Interface script greps energies and gradients to the gradients.dat file in a running directory with the following structure (I states and N atoms):  

```bash
energy ground  
energy 1ex state  
...
...
energy i_state  
grad_x(1at)  grad_y(1at)  grad_z(1at)  
...
...  
grad_x(n_at) grad_y(n_at) grad_z(n_at)  
```

Be carefull: the code expects **gradients** to be extracted, not forces (e.g. Gaussian code). If your ab initio code print forces, you can either grep the negative values of forces to get the gradients (F = - grad(E)) or rename your interface script so that the script name contains one the words forces,gaus or g09 ,separated by dot or minus sign (e.g. forces.sh, forces-gaus.sh, forces-abinitiocode.sh etc), SNAFU will then transform forces into gradients forces.
**NOTE** as the SNAFU is PYTHON based, so the states starts from 0 (ground state). Hence, if an ab initio code denotes ground state as 1, there should be modification in the interface skript (see e.g. molpro examples). 

## TODO:
add diabatization scheme: Le Yu, Phys.Chem.Chem.Phys., 2014, 16, 25883; **doi:10.1039/C4CP03498H**  
