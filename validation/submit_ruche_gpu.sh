#!/bin/bash

#SBATCH --job-name=validation_smilei
#SBATCH --output=%x.o%j 
#SBATCH --time=00:20:00 
#SBATCH --ntasks=1                   # Number of MPI processes (= total number of GPU)
#SBATCH --ntasks-per-node=1          # nombre de tache MPI par noeud (= nombre de GPU par noeud)
#SBATCH --gres=gpu:1
#SBATCH --partition=gpua100 #_test

# Load necessary modules
#source /gpfs/users/prouveurc/env_smilei.sh
source /gpfs/users/prouveurc/myhdf5_env_smilei.sh
module load anaconda3/2022.10/gcc-11.2.0

# Run cuda code
LD_PRELOAD=/gpfs/softs/spack/opt/spack/linux-centos7-haswell/gcc-4.8.5/gcc-11.2.0-mpv3i3uebzvnvz7wxn6ywysd5hftycj3/lib64/libstdc++.so.6.0.29 ./smilei input.py
