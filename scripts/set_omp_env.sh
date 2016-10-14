#!/bin/bash

export OMP_NUM_THREADS=$1
export OMP_SCHEDULE=dynamic
export OMP_PROC_BIND=true

