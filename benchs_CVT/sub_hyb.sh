#!/bin/bash

#@ class            = clallmds+
#@ job_name         = cvt_bench
#@ total_tasks      = 8
#@ node             = 4
#@ as_limit         = 14gb
#@ node_usage       = not_shared
#@ wall_clock_limit = 00:30:00
#@ output           = $(job_name).$(jobid).log
#@ error            = $(job_name).$(jobid).err
#@ job_type         = mpich
#@ environment      = COPY_ALL 
#@ notification     = never
#@ queue


# IntelMPI environment
module load intel/15.0.0 intelmpi/5.0.1 hdf5/1.8.16_intel_intelmpi_mt python/anaconda-2.1.0
unset LD_PRELOAD

export OMP_PROC_BIND=true
export OMP_SCHEDULE=DYNAMIC
export OMP_NUM_THREADS=8

export TASKS_PER_NODE=`echo ${LOADL_PROCESSOR_LIST}|sed -e s/" "/\\n/g|uniq|wc -l`
mpirun -print-rank-map -np ${LOADL_TOTAL_TASKS} -ppn ${TASKS_PER_NODE} ../smilei tst2d_3_plasma_collision.py

