***MCSCF gradients
memory, 1000, m;
file, 2, input.wfu, unknown
PUNCH, input.pun, new
gprint, orbital, civector
orient,noorient;
symmetry,nosym;
Angstrom

geometry=../mini.dat
basis=6-31g

!-we need to get rid of the SAMC records in file 2 (input.wfu,restart file)
!-otherwise, the forces and NACME are wrong for the following reason
!-cpmscf will not rewrite itself it but ather write into following records
!-but the subsequent call to forces would read from old records -> wrong numbers
!-we use file 2 for forces and NACME due to df-casscf

data,truncate,2,3000  !truncate dumpfile after reference (first geometry)

if (lastorb.ne.MCSCF)then
    {hf;wf,9,0,1}
    {multi;
     occ,5;
     closed,3;
     wf,9,0,1;
     state,2;
     maxiter,40;
     orbital,2101.2 !Orbital dumprecord at reference geometry
     save,ci=2501.2}
endif

data, copy, 2101.2, 3000.2

{multi;
occ,5;
closed,3;
WF,9,0,1;
state,2;
maxiter,2;
orbital,2101.2
ciguess,2501.2
save,ci=2501.2
diab,3000.2,save=2101.2,method=-1
}

if (status.lt.0) then
   text, MCSCF failed to converge.
   text, Attempting uncoupled iterations.
   text, Enlarging PSPACE.
   {multi;
   occ,5;
   closed,3;
   WF,9,0,1;
   ! Info about pspace: https://www.molpro.net/info/2015.1/doc/manual/node244.html
   ! uncomment in case of convergence difficulties...
   pspace, 2;
   state,2;
   maxiter,40;
   orbital,2101.2
   ciguess,2501.2
   save,ci=2501.2
   diab,3000.2,save=2101.2,method=-1
   {iterations
   do,uncouple,1,to,5}
   }
   if (status.lt.0) then
      text, Uncoupled iterations did not help, exiting.
      STOP
   endif

endif


TEXT, Calculating gradients. 
multi;
occ,5;
closed,3;
WF,9,0,1;
state,2;
orbital,2101.2
ciguess,2501.2
save,ci=2501.2
diab,3000.2,save=2101.2,method=-1
DM;  ! calculate dipole moments
cpmcscf, grad, 2.1, save=5100.2, accu=1d-;
forces; samc, 5100.2;
pop; density,2101.2,state=1.1
pop; density,2101.2,state=2.1
