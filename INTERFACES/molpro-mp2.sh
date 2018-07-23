#!/bin/bash
source SetEnvironment.sh MOLPRO 2015  # export MOLPROEXE=$molproroot/bin/molpro
cd ABINITIO

##########SNAFU INPUTS###########################################
abinit_geom_file=$1
natoms=$2
state=$3   # for which state
nstate=$4 # total number of state
step=$5
echo "$2 $3 $4 $5" > lal
input=input
nacaccu=9
####################################################################
basis="6-31g"  # for Pople basis sets and df-casscf, fitting DF basis must be specified manually
               # don't use Dunning basis sets, you won't get NACME
nelectrons=9   # total number of electrons
spin=1         # 0 for singlet, 1 for dublet etc.
nocc=5         # occupied orbitals
nclosed=3      # closed orbitals
memory=1000    # molpro memory in MegaWords (1MW = 8 MB)
multi="multi"  # use  "df-casscf" for density fitting version

#if [ -e ../gradients.dat ];then
#rm -f ../gradients.dat
#fi
###################################################################################
cat > $input.com << EOF
memory, $memory,m;
gprint, orbital,civector
symmetry,nosym
Angstrom

geometry=../$abinit_geom_file

basis=$basis

charge=1;
uhf;
ump2;
forces;

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
   cp $input.com.out $input.com.out.error
   if $( grep -q 'NO CONVER' $input.com.out ) ;then 
      echo "ERROR: Could not converge forces!"
      exit 3
   else
      echo "ERROR during execution of MOLPRO. See $input.com.out"
      exit 2
   fi
fi
cp $input.com.out $input.com.out.old
echo "$step" >> grads.dat
################### Extracting energy 
grep "!UMP2 STATE 1.1 Energy " $input.com.out | awk -F "Energy" '{print $2}' | tail -n $nstate > ../gradients.dat
grep "Numerical gradient for UMP2" -A7 $input.com.out | tail -n$natoms  | awk  '{print $2, $3, $4}'>>  ../gradients.dat

if [[ $? -eq 0 ]];then
exit 0
else 
echo "$? Could not extract energies or gradients."
exit 4
fi
