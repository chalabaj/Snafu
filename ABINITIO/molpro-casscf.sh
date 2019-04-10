#!/bin/bash
source SetEnvironment.sh MOLPRO 2015  # export MOLPROEXE=$molproroot/bin/molpro
cd ABINITIO

##########SNAFU INPUTS###########################################
abinit_geom_file=$1
natoms=$2
state=$3     # which electronic state
nstate=$4    # total number of state
step=$5

input=input
nacaccu=9   # forces accuracy
####################################################################
basis="6-31g"  
               # don't use Dunning basis sets, you won't get NACME
nelectrons=9   # total number of electrons
spin=1         # 0 for singlet, 1 for dublet etc.
nocc=6         # occupied orbitals
nclosed=2      # closed orbitals
memory=1000    # molpro memory in MegaWords (1MW = 8 MB)
multi="multi"  # use  "df-casscf" for density fitting version

#if [ -e ../gradients.dat ];then
#rm -f ../gradients.dat
#fi
###################################################################################
cat > $input.com << EOF
memory, $memory,m;
file, 2, $input.wfu,unknown
PUNCH, $input.pun,new
gprint, orbital,civector
symmetry,nosym
Angstrom

geometry=../$abinit_geom_file

basis=$basis

!-for simple CASSCF, you don't need to modify anything below this

!-we need to get rid of the SAMC records in file 2 (input.wfu,restart file)
!-otherwise, the forces and NACME are wrong for the following reason
!-cpmscf will not rewrite itself it but ather write into following records
!-but the subsequent call to forces would read from old records -> wrong numbers
!-we use file 2 for forces and NACME due to df-casscf

data,truncate,2,5101

if (lastorb.ne.MCSCF)then
   {hf;wf,$nelectrons,0,$spin}
endif

$multi;
occ,$nocc;
closed,$nclosed;
WF,$nelectrons,0,$spin;
state,$nstate;
maxiter,40;
ORBITAL,2140.2;
NOEXTRA;

cpmcscf,grad,$state.1,ACCU=1d-$nacaccu,save=5101.2; 
forces;samc,5101.2;

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
   {iterations
   do,uncouple,1,to,5}
   }
endif
EOF

############################----------MOLPRO JOB-------------------------

export TMPDIR=$PWD/scratch

if [[ ! -d $TMPDIR ]]; then
mkdir $TMPDIR
fi

$MOLPROEXE -s --no-xml-output -I $PWD -W $TMPDIR >& $input.com.out <$input.com

# Check whether all is OK.
# If it is some other error, do nothing. It's up to ABIN to decide what to do.
if [[ $? -ne 0 ]];then
    cp $input.com.out $input.com.out.error.$step
    rm $input.pun 
    $MOLPROEXE -s --no-xml-output -I $PWD -W $TMPDIR >& $input.com.out <$input.com
fi 
if [[ $? -ne 0 ]];then
   cp $input.com.out $input.com.out.error.$step.nowf
   if $( grep "NO CONVER" $input.com.out ) ;then 
      echo "ERROR during execution of MOLPRO. See $input.com.out" > MOLPRO_ERROR
      exit 1
   elif $( grep "ERROR DETECTED" $input.com.out ); then
      echo "ERROR during execution of MOLPRO. See $input.com.out" > MOLPRO_ERROR
      exit 2  
   else
      echo "ERROR during execution of MOLPRO. See $input.com.out" > MOLPRO_ERROR
      exit 3
   fi
fi

#echo "$step" >> grads.dat
################### Extracting energy 
grep "MCSCF STATE [[:alnum:]].1 Energy" $input.com.out | awk -F "Energy" '{print $2}' | tail -n $nstate > ../gradients.dat
grep "GRADIENT," $input.pun | awk -F" " '{print $5" "$6" "$7" "}'>> ../gradients.dat #need space at the end for numpy float reading

#cat ../gradients.dat >> grads.dat

if [[ $? -eq 0 ]];then
exit 0
else 
echo "$? Could not extract energies or gradients."
exit 4
fi
