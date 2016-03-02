# environnement 
 module load intel/15.0.0 openmpi  hdf5/1.8.10_intel_openmpi python
  # 
 # execution 
 export OMP_NUM_THREADS=2
 mpirun -bind-to-socket -np 2 /gpfs1l/gpfshome/mds/staff/mfle/GITLAB/smilei/scripts/../validation/workdirs/smilei /gpfs1l/gpfshome/mds/staff/mfle/GITLAB/smilei/scripts/../benchmarks/tst1d_1_clb_explosion.py >smilei_exe.out 
 exit $?  