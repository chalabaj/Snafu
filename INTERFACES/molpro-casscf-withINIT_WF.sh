#!/bin/bash
source SetEnvironment.sh MOLPRO 2015  # export MOLPROEXE=$molproroot/bin/molpro
cd ABINITIO

##########SNAFU INPUTS###########################################
abinit_geom_file=$1
natoms=$2
let state=$3+1   # which electronic state, SNAFU uses 0 for ground state thus the +1
nstate=$4        # total number of state
step=$5

input=input
nacaccu=9   # forces accuracy
####################################################################
basis="6-31g*"  
               # don't use Dunning basis sets, you won't get NACME
nelectrons=40   # total number of electrons
spin=0         # 0 for singlet, 1 for dublet etc.
nocc=22         # occupied orbitals
nclosed=18      # closed orbitals
memory=300    # molpro memory in MegaWords (1MW = 8 MB)
multi="multi"  # use  "df-casscf" for density fitting version
moldenfile="orbitals${step}.molden"
#if [ -e ../gradients.dat ];then
#rm -f ../gradients.dat
#fi
###################################################################################
cat > $input.com << EOF
memory, $memory,m;
file, 2, wf22_18_5st_631gp.wfu,unknown
PUNCH, $input.pun,new
gprint, orbital,civector
symmetry,nosym
Angstrom

geometry=../$abinit_geom_file

basis=$basis

!data,truncate,2,5101

!if (lastorb.ne.MCSCF)then
!   {hf;wf,$nelectrons,0,$spin}
!endif

{$multi;
occ,$nocc;
closed,$nclosed;
WF,$nelectrons,0,$spin;
state,$nstate;
maxiter,40;
START,2140.2,
ORBITAL,2140.2;
NOEXTRA;
cpmcscf,grad,$state.1,ACCU=1d-$nacaccu,save=5101.2}
if ($step.lt.20) then
put,molden,$moldenfile;
end if
forces;samc,5101.2;

EOF

############################----------MOLPRO JOB-------------------------

export TMPDIR=$PWD/scratch

if [[ ! -d $TMPDIR ]]; then
mkdir $TMPDIR
if [[ $step -eq 0 ]]; then
cp wf22_18_5st_631gp.wfu $TMPDIR
fi
fi

$MOLPROEXE -s --no-xml-output -I $PWD -W $TMPDIR >& $input.com.out <$input.com

if [[ $step  -lt 20 ]];then
cp $input.com.out $input.com.out_step${step}
fi

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
