#!/bin/bash
#SBATCH --job-name="install-test"
#SBATCH --nodes=1
#SBATCH --partition=skylake
#SBATCH --ntasks-per-node=32
##


config_dir=configs

spinner_auto_md -np $SLURM_NTASKS total.yaml
spinner_nnp_train total.yaml
configure_csp total.yaml $config_dir
mpirun -np $SLURM_NTASKS spinner_csp $config_dir/final_spinner.yaml
