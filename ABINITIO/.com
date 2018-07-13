***MCSCF gradients
memory, 1000, m;
file, 2, .wfu, unknown
PUNCH, .pun, new
gprint, orbital, civector
Angstrom

geometry=../
basis=6-31g*

!-we need to get rid of the SAMC records in file 2 (input.wfu,restart file)
!-otherwise, the forces and NACME are wrong for the following reason
!-cpmscf will not rewrite itself it but ather write into following records
!-but the subsequent call to forces would read from old records -> wrong numbers
!-we use file 2 for forces and NACME due to df-casscf

data,truncate,2,3000  !truncate dumpfile after reference (first geometry)

if (lastorb.ne.MCSCF)then
   {hf;wf,22,0,0}
   multi;
   occ,14;
   closed,8;
   WF,22,0,0;
   state,;
   maxiter,40;
   orbital,2101.2 !Orbital dumprecord at reference geometry
   save,ci=2501.2
endif

data, copy, 2101.2, 3000.2

{multi;
occ,14;
closed,8;
WF,22,0,0;
state,;
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
   occ,14;
   closed,8;
   WF,22,0,0;
   ! Info about pspace: https://www.molpro.net/info/2015.1/doc/manual/node244.html
   ! uncomment in case of convergence difficulties...
   pspace, 2;
   state,;
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
occ,14;
closed,8;
WF,22,0,0;
state,;
orbital,2101.2
ciguess,2501.2
save,ci=2501.2
diab,3000.2,save=2101.2,method=-1
DM;  ! calculate dipole moments

cpmcscf, grad, .1, save=5100.2, accu=1d-;
forces; samc, 5100.2;
