#!/bin/sh
#PBS -N job_name
#PBS -q batch
#PBS -l nodes=xx:ppn=xx
#PBS -l walltime=480:00:00

cd $PBS_O_WORKDIR
cat $PBS_NODEFILE > nodefile
NPROC=`wc -l < $PBS_NODEFILE`

cd src
mpirun -np $NPROC python3 main.py mg2si1o4.yaml >& ../stdout.x

# if you don't want to print output (mostly useless)
# run such as 
# mpirun -np $NPROC python3 main.py mg2si1o4.yaml > /dev/null 
