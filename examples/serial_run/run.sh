#!/bin/bash
#SBATCH --job-name="spinner"
#SBATCH --partition=csc2          # Partition name (skylake)
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=24         # Cores per node
##

config_dir=configs

spinner_auto_md -np $SLURM_NTASKS total.yaml
spinner_nnp_train total.yaml
configure_csp total.yaml $config_dir
mpirun -np $SLURM_NTASKS spinner_csp $config_dir/final_spinner.yaml
