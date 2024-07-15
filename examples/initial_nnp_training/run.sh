#!/bin/bash
#SBATCH --job-name="spinner-nnp_train"
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --partition=gpu
##


spinner_nnp_train simplenn.yaml
