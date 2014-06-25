#!/bin/bash

if [ ! -d "tests" ]; then
  echo "Execute this script from Smilei root directory, must contain tests directory"
  return
fi

VALIDIR=${PWD}/Validation

mkdir -p ${VALIDIR}/data

cp tests/tst?_*.in                    ${VALIDIR}/data
cp tests/tst2d_0_refl_imm_ions*.in    ${VALIDIR}/data
cp tests/tst2d_1_ion_clb_explosion.in ${VALIDIR}/data

cd ${VALIDIR}
ln -s ../src/smilei
 
for file in $(ls data/tst*.in) ; do
  file_=${VALIDIR}/${file%.in}
  ROOTESTDIR=${file_/data\/tst/TEST} 
  for nprocs in 1 2 ; do
    TESTDIR=${ROOTESTDIR}_${nprocs}_cpu
    mkdir -p ${TESTDIR}
    cd ${TESTDIR}
    ls ../${file}
    # OpenMPI options
    mpirun -np ${nprocs} -bycore -bind-to-core ../smilei ../${file} 2>&1 | tee -a smilei.log
  done
done 

cd ${VALIDIR}/..

