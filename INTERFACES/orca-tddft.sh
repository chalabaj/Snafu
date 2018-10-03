#!/bin/bash
cd ABINITIO
printenv >ll
source SetEnvironment.sh ORCA 4.0.0
export MKL_NUM_THREADS=1
########## SNAFU INPUTS ###########################################
abinit_geom_file=$1
natoms=$2
state=$3   # for which state
let nstate=$4-1 # total number of state minus ground state which is calculated separately
step=$5
input=input.com
if [[ -f $input ]];then
rm $input
fi
########## ORCA PARAMETERS ###########################################
mem=8000   # memory in MB per core
charge=1 
multiplicity=2
basis="def2-svp def2/jk"         # def2-SVP and 6-31+g** are about the same expensive
DFT="BHANDHLYP $basis RIJCOSX"   # TightSCF
#RI-JK approximations for hybrid-DFT methods-: 
#RIJCOSX approximation for Coulomb and HF Exchange is recommended, but RIJK can be better for small molecules
#E(DFT-TDDFT) with RIJCOSX for wat dimer is 0.0001
#E(DFT-TDDFT) with NORI for wat dimer is 0.0000


if [[ $state -eq 0 ]];then   # ground state
 gstask="! ENGRAD $DFT"
 extask="! ENERGY $DFT"
 gradfile="input.engrad"
else                         # excited state
 gstask="! ENERGY $DFT"
 extask="! ENGRAD $DFT"
 iroot="iroot $state"
 gradfile="input_job2.engrad"
fi

if [[ ! -z "${NSLOTS}" ]];then
nprocs=$NSLOTS
cat > $input << EOF
%maxcore $mem          # memory per core
%pal
   nprocs $nprocs           # number of cores
end
EOF
touch nslostsset
fi
########### END OF INPUTS ###########################################

cat >> $input << EOF
$gstask           # for small systems increase accuracy by: Grid5 FinalGrid6 
%scf
maxiter 200
end
! NOPOP
! MiniPrint
! AUTOSTART         # try to read from previous step
* xyz $charge $multiplicity 
EOF
tail -n$natoms ../$abinit_geom_file >> $input
echo '*' >>$input

cat >> $input << EOF
\$new_job
$extask
%scf
maxiter 200
end
! NOPOP
! moread
! MiniPrint

%tddft
  maxdim 5
  nroots $nstate
  maxcore $mem
  $iroot  
  PrintLevel 3
end
* xyz $charge $multiplicity 
EOF
tail -n$natoms ../$abinit_geom_file >> $input
echo '*' >> $input


$ORCAEXE $input > $input.out
################################

if [[ $? -eq 0 ]];then
   cp $input.out $input.out.old
   rm *.tmp
else
   echo "WARNING from ORCA calculation probably failed."
   echo "See $input.out.error_old"
   echo "TRYING TO RESTART JOB WITHOUT PREVIOUS WF" 
   cp $input.out $input.out.error_old
   mv input.com restart.inp
   rm input.*
   mv restart.inp input.com 
   $ORCAEXE $input > $input.out
   if [[ $? -eq 0 ]];then
      cp $input.out $input.out.old
      rm *.tmp
   else
      echo "WARNING from ORCA calculation probably failed."
      echo "See $input.out.error"
      exit 1
      fi
fi


### EXTRACTING ENERGY AND FORCES
# need space at the end for numpy float reading
# rewrite the file

gse=$(grep "E(SCF)" input.com.out | awk {'print $3'})
echo " $gse" > ../gradients.dat  
grep "E\[ [0-9]\]" input.com.out | tail -n$nstate  | awk  -v gsener="$gse" '{sum = $4 + gsener; printf " %.9f \n",sum}' >> ../gradients.dat

awk -v natom="$natoms" '{
if ($4=="gradient") 
    {getline;
	   for(i=1;i<=natom;i++) 
         {getline;x=$1;getline;y=$1;getline;print x, y, $1," "}}}' $gradfile >> ../gradients.dat 

if [[ $? -eq 0 ]];then
exit 0
else 
echo "$? Could not extract energies or gradients."
exit 4
fi
