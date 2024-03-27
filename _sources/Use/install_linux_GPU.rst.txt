
Install on Linux for GPU
-----------------------------

First, make sure you have a recent version of CMAKE, and the other libraries
to compile Smilei on CPU as usual. In particular, for this example, you
need GCC <= 12.

Make a directory to store all the nvidia tools. We call it $NVDIR:

.. code:: bash

  cd $NVDIR

The first step is download ``nvhpc`` and install it:

.. code:: bash

  wget https://developer.download.nvidia.com/hpc-sdk/23.11/nvhpc_2023_2311_Linux_x86_64_cuda_12.3.tar.gz
  tar xpzf nvhpc_2023_2311_Linux_x86_64_cuda_12.3.tar.gz
  nvhpc_2023_2311_Linux_x86_64_cuda_12.3/install # install in $NVDIR when asked
  
  rm nvhpc_2023_2311_Linux_x86_64_cuda_12.3.tar.gz

Then compile ``hdf5``:

.. code:: bash

  export PATH=$NVDIR/Linux_x86_64/23.11/compilers/bin:$PATH
  export PATH=$NVDIR/Linux_x86_64/23.11/comm_libs/mpi/bin:$PATH
  
  wget https://github.com/HDFGroup/hdf5/releases/download/hdf5-1_14_2/hdf5-1_14_2.tar.gz
  tar xpzf hdf5-1_14_2.tar.gz
  
  cd hdfsrc/
  mkdir build
  cd build
  cmake -DCMAKE_C_COMPILER=`which mpicc` -DCMAKE_INSTALL_PREFIX=$NVDIR/hdfsrc/install -DHDF5_ENABLE_PARALLEL=ON ..
  make
  make install

Create a file ``nvidia_env.sh`` containing the following commands:

.. code:: bash

  export BUILD_DIR=build_nvidia
  
  export PATH=$NVDIR/Linux_x86_64/23.11/compilers/bin:$PATH
  export PATH=$NVDIR/Linux_x86_64/23.11/comm_libs/mpi/bin:$PATH
  
  export HDF5_ROOT_DIR=$NVDIR/hdfsrc/install/
  export LD_LIBRARY_PATH=$HDF5_ROOT_DIR/lib
  
  export LDFLAGS="-acc=gpu -gpu=ccnative -cudalib=curand "
  export CXXFLAGS="-acc=gpu -gpu=ccnative,fastmath -std=c++14 -lcurand -Minfo=accel -w -D__GCC_ATOMIC_TEST_AND_SET_TRUEVAL=1 -I$NVDIR/Linux_x86_64/23.11/math_libs/include/"
  export GPU_COMPILER_FLAGS="-O3 --std c++14 -arch=sm_86 --expt-relaxed-constexpr --compiler-bindir $NVDIR -I$NVDIR/Linux_x86_64/23.11/comm_libs/12.3/openmpi4/openmpi-4.1.5/include/ -I$NVDIR/hdfsrc/install/include/"
  
  export SMILEICXX_DEPS=g++
  
  export SLURM_LOCALID=0

To compile Smilei:

.. code:: bash

  source nvidia_env.sh
  make config="gpu_nvidia"

To run:

.. code:: bash

  source nvidia_env.sh
  smilei namelist.py
