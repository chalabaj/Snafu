#!/bin/bash

cd ABINITIO

##########SNAFU INPUTS###########################################
abinit_geom_file=$1
geometry=../$abinit_geom_file

natoms=$2
let natom2=$natoms+2

state=$3   # for which state
nstate=$4 # total number of state
input="input"

source SetEnvironment.sh GAUSSIAN
CHKFILE="$PWD/checkpoint.chk"

#-----------NOTES FOR GAUSSIAN INPUT-------------------
# (guess=(read,Tcheck) reads wavefunction from checkpoint file, if present
# SCF=XQC turns on automatically quadratic convergence (SCF=QC), when we do not converge

# FOR METHODS OTHER THAN HF or DFT, you have to change greping of energies at the end of this script!!!! 

#-----USER SETUP---------------------------------------------------------
cat > $input.com << EOF
\$rungauss
%chk=$CHKFILE
%Mem=500Mb
%Nproc=2
#BLYP/6-31g* Force  guess=(Read,TCheck)

comment

0 1
EOF

########END OF USER MODIFICATIONS###################

cat $geometry >> $input.com
echo ' ' >>$input.com

#######GAUSSIAN STUFF###########
export GAUSS_EXEDIR="$gaussroot"
export GAUSS_SCRDIR="$PWD/scratch"
export scrhome=$GAUSS_SCRDIR
export LD_LIBRARY_PATH=$GAUSS_EXEDIR
mkdir -p $GAUSS_SCRDIR

$GAUSSEXE $input.com

if [[ $? -eq 0 ]];then
   cp $input.log $input.log.old
else
   cp $input.log $input.log.error
   echo "WARNING from r.g09: G09 job probably failed."
   echo " See $input.log.error"
   exit 3
fi

/bin/rm -rf /$GAUSS_SCRDIR
/bin/rm -rf core
################################

### EXTRACTING ENERGY AND FORCES
grep -e 'SCF Done' -e 'EUMP2' $input.log |tail -1|awk '{if ($1=="SCF"){print $5} else if ($1=="E2"){print $6}}' > ../gradients.dat

grep -A$natom2 "Forces (Hartrees/Bohr)" $input.log | tail -n$natoms | awk '{print $3, $4, $5}' >> ../gradients.dat


if [[ $? -eq 0 ]];then
exit 0
else 
echo "$? Could not extract energies or gradients."
exit 4
fi
