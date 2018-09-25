#!/bin/bash
#cd ABINITIO

source SetEnvironment.sh ORCA 4.0.0
##########SNAFU INPUTS###########################################
abinit_geom_file=$1
natoms=$2
state=$3   # for which state
let nstate=$4-1 # total number of state minus ground state which is calculated separately
step=$5
input=input.com
if [[ -f $input ]];then
rm $input
fi
#SYSTEM:
mem=8000   # memory in MB per core
charge=0 
multiplicity=1
basis=LANL2DZ         # def2-SVP and 6-31+g** are about the same expensive
DFT="B3LYP RIJCOSX"   # TightSCF 
#RI-JK approximations for hybrid-DFT methods-: 
#RIJCOSX approximation for Coulomb and HF Exchange is recommended, but RIJK can be better for small molecules
#E(DFT-TDDFT) with RIJCOSX for wat dimer is 0.0001
#E(DFT-TDDFT) with NORI for wat dimer is 0.0000


if [[ $state -eq 0 ]];then
 gstask="! ENGRAD $DFT"
 extask="! ENERGY $DFT"
else
 gstask="! ENERGY $DFT"
 extask="! ENGRAD $DFT"
 iroot="iroot $state"
fi

if [[ -z {$NSLOTS} ]];then
nprocs=$NSLOTS
cat > $input << EOF
%maxcore $mem          # memory per core
%pal
   nprocs $nprocs           # number of cores
end
EOF
fi
###########END OF INPUTS

cat >> $input << EOF
%basis
Basis "$basis"        # The orbital expansion basis set
end

$gstask           # for small systems increase accuracy by: Grid5 FinalGrid6 
! AUTOSTART         # try to read from previous step
* xyz $charge $multiplicity 
EOF
tail -n$natoms ../$abinit_geom_file >> $input
echo '*' >>$input

cat >> $input << EOF
\$new_job
$extask
! moread
%basis
 Basis "$basis"        # The orbital expansion basis set
end
%tddft
  maxdim 5
  nroots $nstate
  maxcore $mem
  $iroot  
end
* xyz $charge $multiplicity 
EOF
tail -n$natoms ../$abinit_geom_file >> $input
echo '*' >> $input


$ORCAEXE $input > $input.out
################################
if [[ $? -eq 0 ]];then
   cp $input.out $input.out.old
else
   echo "WARNING from r.orca: ORCA calculation probably failed."
   echo "See $input.out.error" 
   cp $input.out $input.out.error
fi

### EXTRACTING ENERGY AND FORCES
awk -v natom="$natoms" '{
if ($4=="gradient") {getline;
	for(i=1;i<=natom;i++) {
		getline; x=$1;getline;y=$1;getline;print x,y,$1
 }}}' *.engrad > ../engrad.dat 

