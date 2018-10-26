
Install dependencies on Linux
-----------------------------


Fedora
^^^^^^^^^

.. code-block:: bash

  dnf install -y gcc-c++ hdf5-openmpi hdf5-openmpi-devel openmpi-devel git which findutils python python-devel
  dnf install -y h5py ipython python2-pint sphinx python2-matplotlib



Debian (Ubuntu, Mint etc...)
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

1. Install these packages

  .. code-block:: bash
  
    sudo apt-get install python-h5py ipython python-pint python-sphinx python-matplotlib python-dev  python-numpy

2. Since the system ``openmpi`` is not compiled with
   ``--enable-mpi-thread-multiple``, a manual installation is required.
   Add the following lines to your `~/.bashrc` or `~/.bash_profile` file
   (You may choose any ``${INSTALL_DIR}``)

  .. code-block:: bash

    export INSTALL_DIR=/usr/local
    export PATH=${INSTALL_DIR}/openmpi/bin:${PATH}
    export LD_LIBRARY_PATH=${INSTALL_DIR}/openmpi/lib:${LD_LIBRARY_PATH}
    export PATH=${INSTALL_DIR}/hdf5/bin:${PATH}
    export LD_LIBRARY_PATH=${INSTALL_DIR}/hdf5/lib:${LD_LIBRARY_PATH}
    export HDF5_ROOT_DIR=${INSTALL_DIR}/hdf5

3. Restart your terminal

4. Download `OpenMPI <https://www.open-mpi.org/software/ompi>`_ and install.

  .. code-block:: bash
  
    tar zxvf openmpi-*.*.*.tar.gz
    cd openmpi-*.*.*
    ./configure --prefix=${INSTALL_DIR}/openmpi --enable-mpi-thread-multiple --enable-mpirun-prefix-by-default
    make
    sudo make install

5. Restart your terminal

6. Download `HDF5 <https://portal.hdfgroup.org/display/support/Downloads>`_ and install

  .. code-block:: bash

    tar zxvf hdf5-*.*.*.tar.gz
    cd hdf5-*.*.*
    ./configure --prefix=${INSTALL_DIR}/hdf5 --enable-parallel --with-pic --enable-linux-lfs --enable-shared --enable-build-mode=production --disable-sharedlib-rpath --enable-static CC=mpicc FC=mpif90
    make
    sudo make install

