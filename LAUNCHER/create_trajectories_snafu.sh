#!/bin/bash
# Script rewritten from ABINs create_trajetories.sh
# This script generates and executes a set of dynamical trajectories using SNAFU.
# Initial geometries are taken sequentially from a XYZ movie file.
# Initial velocities are optiononaly taken sequentially from a XYZ file.
# The trajectories are executed and stored in $folder.
#---------------------------------------------------------------------------------

#######-----SETUP---#############
movie=movpick.xyz      # PATH TO a XYZ movie with initial geometries
veloc=velpick.xyz      # PATH to XYZ initial velocities, leave blank if you do not have them
isample=1	             # initial number of traj
nsample=100	           # number of trajectories
folder=SNAFU           # Name of the folder with trajectories
inputdir=TEMPLATE-$folder   # Directory with input files for ABIN
abin_input=$inputdir/input.in   # main input file for ABIN
launch_script=$inputdir/SNAFUS	# this is the file that is submitted by qsub - copy of file SNAFU/SNAFUS
submit="qsub -q nq-16-8 -V -cwd  " # comment this line if you don't want to submit to queue yet
rewrite=0                          # if =1 -> rewrite trajectories that already exist
jobs=8                             # number of batch jobs to submit. Trajectories will be distributed accordingly.

# Number of atoms is determined automatically from input.in
natom=$(awk -F"[! ,=]+" '{if($1=="natoms")print $2}' $abin_input) #number of atoms
molname=$folder      # Name of the job in the queue
##########END OF SETUP##########


function Folder_not_found {
   echo "Error: Folder $1 does not exists!"
   exit 1
}

function File_not_found {
   echo "Error: File $1 does not exists!"
   exit 1
}

function Error {
   echo "Error from command $1. Exiting!"
   exit 1
}

if [[ ! -d "$inputdir" ]];then
   Folder_not_found "$inputdir"
fi

if [[ ! -f "$movie" ]];then
   File_not_found "$movie"
fi

if [[ ! -z "$veloc" ]] && [[ ! -f "$veloc" ]];then
   File_not_found "$veloc"
fi

if [[ ! -e $abin_input ]];then
   File_not_found "$abin_input"
fi

if [[ ! -e "$launch_script" ]];then
   File_not_found "$launch_script"
fi

if [[ -e "geom.in" ]] || [[ -e "restart.xyz" ]];then
   echo "Error: Files geom.in and/or restart.xyz were found here."
   echo "Please remove them."
   exit 1
fi

#   -------------------------------------------------------------------------------------

echo "Number of atoms = $natom"

let natom2=natom+2
lines=$(cat $movie | wc -l)
geoms=$(expr $lines / $natom2)
if [[ $nsample -gt $geoms ]];then
   echo "ERROR: Number of geometries ($geoms) is smaller than number of samples($nsample)."
   echo "Change parameter \"nsample\"."
   exit 1
fi


# determine number simulations per job
let nsimul=nsample-isample+1
if [[ $nsimul -le $jobs ]];then
   remainder=0
   injob=1
   jobs=$nsimul
else
   let injob=nsimul/jobs  #number of simulations per job
   # determine the remainder and distribute it evenly between jobs
   let remainder=nsimul-injob*jobs
fi


j=1
i=$isample

#--------------------------------------------------------------------------------

mkdir -p $folder
cp "$abin_input" $folder

let offset=natom2*isample-natom2

while [[ $i -le "$nsample" ]];do

   let offset=offset+natom2     

   if [[ -d "$folder/TRAJ.$i" ]];then
      if [[ "$rewrite" -eq "1" ]];then

         rm -r $folder/TRAJ.$i ; mkdir $folder/TRAJ.$i
         rm -f $folder/$molname.$isample.$j.sh

      else

         echo "Trajectory number $i already exists!"
         echo "Exiting..."
         exit 1

      fi

   else

      mkdir $folder/TRAJ.$i

   fi

   cp -r $inputdir/* $folder/TRAJ.$i


#--Now prepare geom.in (and possibly veloc.in)

   head -$offset $movie | tail -$natom2 > geom
   mv geom $folder/TRAJ.$i/geom.in 
   
   if [[ ! -z "$veloc" ]];then
      head -$offset $veloc | tail -$natom2 > veloc.in
      mv veloc.in $folder/TRAJ.$i/
   fi


## Now prepare launch scripts 

   echo "cd TRAJ.$i" >> $folder/$molname.$isample.$j.sh
   echo "./SNAFUS" >> $folder/$molname.$isample.$j.sh                   
   echo "cd $PWD/$folder" >> $folder/$molname.$isample.$j.sh

#--Distribute calculations evenly between jobs for queue
   if [[ $remainder -le 0 ]];then
      let ncalc=injob
   else
      let ncalc=injob+1 
   fi
   if [[ `expr \( $i - $first + 1 \) % $ncalc` -eq 0 ]] && [[ $j -lt $jobs ]]; then
      let j++
      let remainder--
   fi
#---------------------------------------------------------------------------

   let i++

done

# Submit jobs

if [[ ! -z "$submit" ]];then
   cd $folder
   while [[ $k -le $j ]]
   do
      if [[ -f $molname.$isample.$k.sh ]];then
         $submit -V -cwd $molname.$isample.$k.sh
      fi
      let k++
   done
fi

