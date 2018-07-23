memory, 1000,m;
file, 2, input.wfu,unknown
PUNCH, input.pun,new
gprint, orbital,civector
symmetry,nosym
Angstrom

geometry=../abinit_geom.xyz

basis=6-31g

!-for simple CASSCF, you don't need to modify anything below this

!-we need to get rid of the SAMC records in file 2 (input.wfu,restart file)
!-otherwise, the forces and NACME are wrong for the following reason
!-cpmscf will not rewrite itself it but ather write into following records
!-but the subsequent call to forces would read from old records -> wrong numbers
!-we use file 2 for forces and NACME due to df-casscf

data,truncate,2,5101

if (lastorb.ne.MCSCF)then
   {hf;wf,9,0,1}
endif

multi;
occ,5;
closed,3;
WF,9,0,1;
state,2;
maxiter,40;
ORBITAL,2140.2;
NOEXTRA;

cpmcscf,grad,1.1,ACCU=1d-9,save=5101.2; 
forces;samc,5101.2;

