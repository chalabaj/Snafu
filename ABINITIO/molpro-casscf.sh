#!/bin/bash
source SetEnvironment.sh MOLPRO 2015   # export MOLPROEXE=$molproroot/bin/molpro
cd ABINITIO

##########SNAFU INPUTS###########################################
abinit_geom_file=$1
natoms=$2
state=$3   # for which state
nstate=$4 # total number of state
input=input
####################################################################
basis="6-31g"  # for Pople basis sets and df-casscf, fitting DF basis must be specified manually
               # don't use Dunning basis sets, you won't get NACME
nelectrons=9   # total number of electrons
spin=1         # 0 for singlet, 1 for dublet etc.
nocc=5         # occupied orbitals
nclosed=3      # closed orbitals
memory=1000    # molpro memory in MegaWords (1MW = 8 MB)
multi="multi"  # use  "df-casscf" for density fitting version

if [ -e ../gradients.dat ];then
rm -f ../gradients.dat
fi

# MOLPRO-CASSCF
# wavefunction passed between steps via input.wfu
# 3rd line is needed, we take forces from input.pun
cat > $input.com << EOF
***MCSCF gradients
memory, $memory, m;
file, 2, $input.wfu, unknown
PUNCH, $input.pun, new
gprint, orbital, civector
!orient,noorient;
symmetry,nosym;
Angstrom

geometry=../$abinit_geom_file
basis=$basis

!-we need to get rid of the SAMC records in file 2 (input.wfu,restart file)
!-otherwise, the forces and NACME are wrong for the following reason
!-cpmscf will not rewrite itself it but ather write into following records
!-but the subsequent call to forces would read from old records -> wrong numbers
!-we use file 2 for forces and NACME due to df-casscf

data,truncate,2,5100

if (lastorb.ne.MCSCF)then
    {hf;wf,$nelectrons,0,$spin}

endif

{$multi;
occ,$nocc;
closed,$nclosed;
WF,$nelectrons,0,$spin;
state,$nstate;
maxiter,40;
DM;  ! calculate dipole moments
orbital,2101.2;
ciguess,2501.2
save,ci=2501.2
}

if (status.lt.0) then
   text, MCSCF failed to converge.
   text, Attempting uncoupled iterations.
   text, Enlarging PSPACE.
   {$multi;
   occ,$nocc;
   closed,$nclosed;
   WF,$nelectrons,0,$spin;
   ! Info about pspace: https://www.molpro.net/info/2015.1/doc/manual/node244.html
   ! uncomment in case of convergence difficulties...
   pspace, 2;
   state,$nstate;
   maxiter,40;
   ciguess,2501.2
   save,ci=2501.2
   {iterations
   do,uncouple,1,to,5}
   }
   if (status.lt.0) then
      text, Uncoupled iterations did not help, exiting.
      STOP
   endif

endif

TEXT, Calculating gradients. 
$multi;
occ,$nocc;
closed,$nclosed;
WF,$nelectrons,0,$spin;
state,$nstate;
ciguess,2501.2
save,ci=2501.2
DM;  ! calculate dipole moments
EOF

rec=5100
# Store corresponding cpmscf and forces commands
cpmscf="cpmcscf, grad, $state.1, save=$rec.2, accu=1d-$accu;"
echo "${cpmscf}" >> $input.com
forces="forces; samc, $rec.2;"
echo ${forces} >> $input.com

# Calculate Mullican charges
for ((ist1=1;ist1<=nstate;ist1++)) {
echo "pop; density,2140.2,state=$ist1.1" >> $input.com
}

#----------MOLPRO JOB-------------------------
export TMPDIR=$PWD/scratch

if [[ ! -d $TMPDIR ]]; then
mkdir $TMPDIR
fi

$MOLPROEXE -s --no-xml-output -I $PWD -W $TMPDIR >& $input.com.out <$input.com

cp $input.com.out $input.com.out.old
# Check whether all is OK.
# If it is some other error, do nothing. It's up to ABIN to decide what to do.
if [[ $? -ne 0 ]];then

   cp $input.com.out $input.com.out.error
   if $( grep -q 'NO CONVER' $input.com.out ) ;then 

      echo "ERROR: Could not converge forces!"
      exit 3

   else

      echo "ERROR during execution of MOLPRO. See $input.com.out"
      exit 2

   fi
fi

#####################################################################


# Extracting energy 
grep "MCSCF STATE [[:alnum:]].1 Energy" $input.com.out | awk -F "Energy" '{print $2}' | tail -n $nstate > ../gradients.dat
#grep "MCSCF STATE $state.1 Energy" $input.com.out | awk -F "Energy" '{print $2}' | tail -n $nstate >> lal

# Extracting GRADIENT
# Should work for both CASPT2 and CASSCF gradients
grep "GRADIENT," $input.pun | awk -F" " '{print $5,$6,$7}'>> ../gradients.dat

if [[ $? -eq 0 ]];then
exit 0
else 
echo "$? Could not extract energies or gradients."
exit 4
fi
