#!/bin/bash
cd ABINITIO
abinit_geom_file=$1
natoms=$2
state=$3
nstates=$4

echo "la" > abinit.out

grep "la" out
if [ "$?" == "0" ]; then
 exit 0
else
 exit 6
fi

