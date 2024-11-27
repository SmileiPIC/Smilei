#!/bin/bash

#SBATCH --job-name=compilation_smilei
#SBATCH --output=%x.o%j 
#SBATCH --time=00:20:00 
#SBATCH --ntasks=40                   # Number of MPI processes (= total number of GPU)
#SBATCH --ntasks-per-node=40          # nombre de tache MPI par noeud (= nombre de GPU par noeud), max 80 in theory
#SBATCH --partition=cpu_short

# Load necessary modules
#source /gpfs/users/prouveurc/env_smilei.sh
source /gpfs/users/prouveurc/myhdf5_env_smilei.sh
module load anaconda3/2022.10/gcc-11.2.0

# Run cuda code
make -j 40  machine="ruche_gpu2" config="gpu_nvidia noopenmp verbose" 
