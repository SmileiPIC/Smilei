
Install dependencies on Linux
-----------------------------

Here we present the packages you need to install in order to be able to compile Smilei. Please be aware that distribution change quite often the package names. As 
result this guide could partially or totally outdated. In case you find some error, please fill a `github issue <https://github.com/SmileiPIC/Smilei/issues/new?assignees=&labels=installation&template=installation-errors.md&title=>`_ . 


ArchLinux
^^^^^^^^^

.. code-block:: bash

	sudo pacman -S git hdf5-openmpi python-numpy python-sphinx python-h5py-openmpi python-matplotlib python-pint make gcc


Fedora
^^^^^^

.. code-block:: bash

	sudo dnf install gcc-c++ git hdf5-openmpi hdf5-openmpi-devel openmpi-devel python python-devel python3-h5py ipython python3-pint python3-sphinx python3-matplotlib

Add the following lines to your ``~/.bashrc`` or ``~/.bash_profile`` file

.. code-block:: bash

	module load mpi


Debian or Ubuntu
^^^^^^^^^^^^^^^^

.. code-block:: bash

	sudo apt-get install git python3-h5py python3-ipython python3-pint python3-sphinx python3-matplotlib python3-dev python3-numpy build-essential gcc libhdf5-openmpi-dev

Add the following lines to your ``~/.bashrc`` or ``~/.bash_profile`` file

.. code-block:: bash

	export PYTHONEXE=python3 
	export HDF5_ROOT_DIR=/usr/lib/x86_64-linux-gnu/hdf5/openmpi 


Troubleshooting:
^^^^^^^^^^^^^^^^

Besides Python Smilei need a quite recent mpi (with mpi-thread-multiple enabled) and a parallel hdf5 library. 
In case you system doe not provide them, here is a (non exhaustive) help to instlal them:


1. If your system ``openmpi`` is not compiled with ``--enable-mpi-thread-multiple``, a manual installation is required.
   Add the following lines to your ``~/.bashrc`` or ``~/.bash_profile`` file
   (You may choose any ``${INSTALL_DIR}``)

  .. code-block:: bash

    export INSTALL_DIR=/usr/local
    export PATH=${INSTALL_DIR}/openmpi/bin:${PATH}
    export LD_LIBRARY_PATH=${INSTALL_DIR}/openmpi/lib:${LD_LIBRARY_PATH}
    export PATH=${INSTALL_DIR}/hdf5/bin:${PATH}
    export LD_LIBRARY_PATH=${INSTALL_DIR}/hdf5/lib:${LD_LIBRARY_PATH}
    export HDF5_ROOT_DIR=${INSTALL_DIR}/hdf5

2. Restart your terminal

3. Download `OpenMPI <https://www.open-mpi.org/software/ompi>`_ and install.

  .. code-block:: bash
  
    tar zxvf openmpi-*.*.*.tar.gz
    cd openmpi-*.*.*
    ./configure --prefix=${INSTALL_DIR}/openmpi --enable-mpi-thread-multiple --enable-mpirun-prefix-by-default
    make
    sudo make install

4. Restart your terminal

5. Download `HDF5 <https://portal.hdfgroup.org/display/support/Downloads>`_ and install

  .. code-block:: bash

    tar zxvf hdf5-*.*.*.tar.gz
    cd hdf5-*.*.*
    ./configure --prefix=${INSTALL_DIR}/hdf5 --enable-parallel --with-pic --enable-linux-lfs --enable-shared --enable-build-mode=production --disable-sharedlib-rpath --enable-static CC=mpicc FC=mpif90
    make
    sudo make install

