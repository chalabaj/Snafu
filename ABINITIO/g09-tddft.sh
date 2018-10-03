#!/bin/bash
source SetEnvironment.sh GAUSSIAN
#cd ABINITIO
########## USER INPUT ###########################
DFT=BLYP              
charge=0 
multiplicity=1
basis="6-31g*"         # def2-SVP and 6-31+g** are about the same expensive


##########SNAFU INPUTS###########################
abinit_geom_file=$1
geometry=../$abinit_geom_file
natoms=$2
let natom2=$natoms+2
state=$3   # for which state
let nstate=$4-1 # total number of state minus ground state which is calculated separately
input="input"
CHKFILE="checkpoint.chk"

if [[ ! -z "${NSLOTS}" ]];then
let mem=${NSLOTS}*4
else
NSLOTS=1
fi

if [[ $state -eq 0 ]];then   # ground state
 gstask="#P $DFT/$basis* SCF=XQC Force guess=(Read,TCheck)"
 extask="#P $DFT/$basis* SCF=XQC TD=(nstates=$nstate) Geom=Check Guess=Read"
else                         # excited state
 gstask="#P $DFT/$basis* SCF=XQC guess=(Read,TCheck)"
 extask="#P $DFT/$basis* SCF=XQC FORCE TD=(nstates=$nstate,root=$state) Geom=Check Guess=Read"
fi

cat > $input.com << EOF
\$rungauss
%chk=$CHKFILE
%Mem=${mem}Mb
%Nproc=${NSLOTS} 
$gstask

comment

$charge $multiplicity
EOF

########END OF USER MODIFICATIONS###################

cat $geometry >> $input.com
echo ' ' >>$input.com

cat >> $input.com << EOF
--Link1--
%chk=$CHKFILE
%Mem=${mem}Mb
%Nproc=${NSLOTS} 
$extask

comment

$charge $multiplicity

EOF
#######GAUSSIAN STUFF###########
export GAUSS_EXEDIR="$gaussroot"
export GAUSS_SCRDIR="$PWD/scratch"
export scrhome=$GAUSS_SCRDIR
export LD_LIBRARY_PATH=$GAUSS_EXEDIR
if [[ ! -d $GAUSS_SCRDIR ]];then
mkdir -p $GAUSS_SCRDIR
fi
$GAUSSEXE $input.com

if [[ ! $? -eq 0 ]];then
   echo "WARNING from GAUSSIAN calculation probably failed."
   echo "See $input.out_error_old"
   echo "TRYING TO RESTART JOB WITHOUT PREVIOUS WF" 
   cp $input.log err_$input.out_old
   mv $input.com restart.inp
   rm input.*
   mv restart.inp $input.com
   $GAUSSEXE $input.com

   if [[ !$? -eq 0 ]];then
       echo "WARNING from GAUSSIAN calculation probably failed."
       echo "See $input.out.error"
       exit 1
   fi
fi
cp $input.log $input.log_old

/bin/rm -rf /$GAUSS_SCRDIR
/bin/rm -rf core
################################

### EXTRACTING ENERGY AND FORCES
gse1=$(grep "SCF Done:  E([A-Za-z]" $input.log | tail -1 | awk '{print $5}')
gse2=$(grep "SCF Done:  E([A-Za-z]" $input.log | head -1 | awk '{print $5}')
if [[ $gse1 == $gse2 ]];then
    echo "$gse1" > ../gradients.dat  
else
    echo "SCF energies between DFT and TD-DFT does not match."
    exit 1
fi
grep "Convergence achieved on expansion vectors." $input.log -B$nstate | head -n$nstate  | awk  -v gsener="$gse1" '{sum = ($4 / 27.211396) + gsener; printf " %.9f \n",sum}' >> ../gradients.dat

grep -A$natom2 "Forces (Hartrees/Bohr)" $input.log | tail -n$natoms | awk '{print $3, $4, $5}' >> ../gradients.dat




if [[ $? -eq 0 ]];then
exit 0
else 
echo "$? Could not extract energies or gradients."
exit 4
fi
