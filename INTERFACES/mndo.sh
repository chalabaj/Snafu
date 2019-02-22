#!/bin/bash
#source SetEnvironment.sh MNDO 7  # export MOLPROEXE=$molproroot/bin/molpro
MNDOEXE=MNDO
cd ABINITIO

# SNAFU INPUTS #######################################################################
abinit_geom_file=$1
natoms=$2
let state=$3+1     # which electronic state, SNAFU uses 0 for ground state thus the +1
nstate=$4    # total number of state
step=$5
# USER INPUT #########################################################################
method="OM3"
charge=1       # molecular charge
multi=0        # multiplicity
	       # 0 Closed-shell singlet
	       # 1 Open-shell singlet with two singly occupied orbitals, this usually corresponds to an excited singlet state
	       # 2 Doublet
	       # 3 Triplet
naocc=6         # number of active occupied orbitals
naunocc=2       # number of active unoccupied orbitals
nref=1          # number of reference occupations (see also refdef=3 option for automatic nref extension)
		# 0 None. Full CI in the active space.
		# n Chosen number, maximum 3 (usually, can be more with manual selection of occupations)
refdef=3	# definition of reference occupations
		# 0 and nciref=1: SCF configuration
                #       nciref=2: SCF configuration and doubly excited HOMO-LUMO
                #       nciref=3: SCF configuration, singly and doubly excited HOMO-LUMO
		# 3 Starting as mciref=0, then adds further references so that their fraction is at least 85 %
		#   It efficiently adds more references then defined by nref
nexc=2		# maximum excitation level, ignored for nref=0
                # 1 CIS, only single excitations.
                # 2 CISD, up to double excitations.
                # 3 CISDT, up to triple excitations.
                # 4 CISDTQ, up to quadruple excitations.
                # n Up to n-fold excitations
#memory=1000    # memory in

###################################################################################
input="input"
cat > $input.com << EOF
$method jop=-2 igeom=1 iform=1 nsav15=3 kharge=$charge imult=$multi +
kci=5 ici1=$naocc ici2=$naunocc iroot=$nstate lroot=$state ipubo=1 +
nciref=$nref mciref=$refdef levexc=$nexc ioutci=1 kitscf=500 ktrial=11
Step $step, state $state

EOF

awk 'BEGIN{
   elems["H"]=1
   elems["He"]=2
   elems["Li"]=3
   elems["Be"]=4
   elems["B"]=5
   elems["C"]=6
   elems["N"]=7
   elems["O"]=8
   elems["F"]=9
   elems["Ne"]=10
   elems["Na"]=11
   elems["Mg"]=12
   elems["Al"]=13
   elems["Si"]=14
   elems["P"]=15
   elems["S"]=16
   elems["Cl"]=17
   elems["Ar"]=18
}
#conversion of xyz input to dftb geom
{
   print elems[$1], $2, 0, $3, 0, $4, 0
}'  ../$abinit_geom_file >> $input.com

# JOB #####################################################################################

export TMPDIR=$PWD/scratch

if [[ ! -d $TMPDIR ]]; then
mkdir $TMPDIR
fi

$MNDOEXE $input.com

if [[ $? -ne 0 ]] || ! $( grep -q "COMPUTATION TIME" $input.com.out ) || ! $( grep -q "SCF TOTAL ENERGY" $input.com.out ) ;then
   cp $input.com.out $input.com.out.error.$step
   echo "ERROR during execution of MNDO. See $input.com.out" > ERROR
   exit 1
fi

# EXTRACTION ##############################################################################

natom=$(grep "CARTESIAN COORDINATES: NUMAT" fort.15 | awk '{print $5}')
grep "State" $input.com.out | awk '{print $9/27.211396}' > ../gradients.dat
grep -A $natom "CARTESIAN GRADIENT" fort.15 | tail -n +2 | awk 'BEGIN {toAU=1.5936e-3/1.8897259886}{print -$3*toAU, -$4*toAU, -$5*toAU}' >> ../gradients.dat
