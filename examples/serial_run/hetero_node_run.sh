#!/bin/bash
#SBATCH --job-name="spinner-hetero-node"
#SBATCH --partition=csc2 --ntasks-per-node=32
#SBATCH hetjob
#SBATCH --partition=gpu2 --ntasks-per-node=1
##

config_dir=configs

spinner_auto_md -np $SLURM_NTASKS total.yaml
srun --het-group=1 spinner_nnp_train total.yaml
configure_csp total.yaml $config_dir
mpirun -np $SLURM_NTASKS spinner_csp $config_dir/final_spinner.yaml
