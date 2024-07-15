#!/bin/bash
#SBATCH --job-name="spinner-csp"
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=32
#SBATCH --partition=csc2
##


mpirun -np $SLURM_NTASKS spinner_csp final_spinner.yaml
