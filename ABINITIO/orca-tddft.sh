#!/bin/bash
cd ABINITIO

source ../SetEnvironment.sh ORCA 4.0.0
##########SNAFU INPUTS###########################################
abinit_geom_file=$1
natoms=$2
state=$3   # for which state
nstate=$4 # total number of state
step=$5
input=input.com
#SYSTEM:
charge=-3 
multiplicity=6
basis=LANL2DZ
DFTfunc="B3LYP"

mem=8000   # memory in MB per core

if [[ ! -z  {$NSLOTS}]]; then
nprocs=$NSLOTS
else
nprocs=1
fi
###########END OF INPUTS

cat > $input << EOF

%maxcore $mem          # memory per core
%pal 
   nprocs $nprocs           # number of cores
end

%basis
Basis "$basis"        # The orbital expansion basis set
end

! ENGRAD $DFTfunc  RIJCOSX # for small systems increase accuracy by: Grid5 FinalGrid6 
! AUTOSTART TightSCF   # try to read from previous step
* xyz $charge $multiplicity 
EOF
tail -n$natoms ../$abinit_geom_file >> $input
echo '*' >>$input

if [[ $state -eq 0 ]];then

cat >> $input << EOF
$new_job
! $DFTfunc RIJCOSX TightSCF
! moread
%basis
 Basis "LANL2DZ"        # The orbital expansion basis set
end
%tddft
  maxdim 5
  nroots 5
  maxcore 32000
end
* xyz $charge $multiplicity 
EOF
tail -n$natoms ../$abinit_geom_file >> $input
else 
%tddft
  maxdim 5
  nroots $nstate
  iroot $state
  maxcore 32000
end
fi

export TMPDIR=$PWD/scratch

if [[ ! -d $TMPDIR ]]; then
mkdir $TMPDIR
fi

$ORCAEXE $input &> $input.out
################################
if [[ $? -eq 0 ]];then
   cp $input.out $input.out.old
else
   echo "WARNING from r.orca: ORCA calculation probably failed."
   echo "See $input.out.error" 
   cp $input.out $input.out.error
fi

### EXTRACTING ENERGY AND FORCES
awk -v natom="$natom" '{if ($5=="energy") {getline;getline;print $1}
if ($4=="gradient") {getline;
	for(i=1;i<=natom;i++) {
		getline; x=$1;getline;y=$1;getline;print x,y,$1
 }}}' $input.engrad > ../engrad.dat.$ibead

