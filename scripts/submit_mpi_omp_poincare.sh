#!/bin/bash
#@ class            = clallmds+
#@ as_limit         = 14gb
#@ node_usage       = not_shared 
#@ job_name         = long
#@ total_tasks      = 80 
#@ node             = 40      
#@ wall_clock_limit = 20:00:00
#@ output           = $(job_name).$(jobid).log
#@ error            = $(job_name).$(jobid).err
#@ job_type         = mpich
#@ environment      = COPY_ALL
#@ queue

module load intel openmpi
export OMP_NUM_THREADS=8
export OMP_SCHEDULE=dynamic
export OMP_PROC_BIND=true

mpirun -bysocket -bind-to-socket ../src/smilei test_reconnection.in


