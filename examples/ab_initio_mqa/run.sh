#!/bin/bash
#SBATCH --job-name="spinner-abinitio-mqa"
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=32
#SBATCH --partition=csc2
##


spinner_auto_md -np $SLURM_NTASKS auto_md.yaml
