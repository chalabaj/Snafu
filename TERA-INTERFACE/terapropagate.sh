
# Script for launching terachem server and python server for fast calculation of bunch of jobs
# Stepan Srsen
# -V
export PATH="/home/srsen/bin/anaconda3/bin:$PATH"
export LD_LIBRARY_PATH=
source SetEnvironment.sh TERACHEM 1.9-dev
MPIRUN_TERA="$MPIRUN -np 1 "
source SetEnvironment.sh ABIN mpi
MPITYPE=2
port=teraport$$
#rm RecalcGeometriesTERA.sh.*
rm -rf tera2.out* scr

function ifkill {
	if `ps|grep -q $1` ;then
        	kill -9 $1
	fi

}

trap $(ifkill $terapid) SIGUSR2 
trap $(ifkill $terapid) EXIT 

$MPIRUN_TERA $TERAEXE --inputfile=tera2.inp --UseMPI=$MPITYPE --MPIPort=$port >> tera2.out 2>&1 &
terapid=$!
#sleep 15
for i in {1..100}
do
	if grep -q port_name: tera2.out; then
		break
	elif [ $i -eq 100 ]; then
		ifkill $terapid
		echo "The port for the terachem connection was not published in the time limit."
		exit
	fi
	sleep 1
done
port_2=`grep port_name: tera2.out | awk '{print $6}' | tail -1`
export MPI_TERA_PORT=${port_2}
$MPIRUN_TERA python -u tera-propagate.py

if [[ $? -ne 0 ]]; then
    echo "python failed"
    ifkill $terapid
fi

pythonpid=$!
sleep 1
ifkill $terapid
ifkill $pythonpid