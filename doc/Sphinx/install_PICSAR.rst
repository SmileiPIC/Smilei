
The optional Picsar library
-----------------------------

A Pseudo-Spectral Analytical Time Domain solver (PSATD) for Maxwell equations
is available for use in Smilei via the `Picsar <http://picsar.net>`_ library.
The PSATD solver is an FFT-based high-order Maxwell solver.
Picsar uses the `FFTW <http://www.fftw.org>`_ library.

WARNING: This feature is experimental.

Install FFTW
^^^^^^^^^^^^^^^^^^^^

Download and install the `latest version <http://www.fftw.org/>`_ of FFTW

.. code-block:: bash

  tar zxvf fftw-3.3.7.tar.gz
  cd fftw-3.3.7
  configure --prefix INSTALL_DIR --enable-shared --enable-threads --with-openmp --enable-mpi
  make 
  make install

Set a few environment variables, typically in your `~/.bashrc` or `~/.bash_profile` file.

.. code-block:: bash

  export FFTW_LIB_DIR=${INSTALL_DIR}/lib
  export FFTW_INC_DIR=${INSTALL_DIR}/include
  export LD_LIBRARY_PATH=${INSTALL_DIR}/lib:${LD_LIBRARY_PATH} 


Install PICSAR as a library
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

1. Download the `latest version of Picsar <https://bitbucket.org/berkeleylab/picsar.git>`_

  .. code-block:: bash  

    git clone git@bitbucket.org:berkeleylab/picsar.git
    cd picsar/

2. Set library flags to compile picsar as a library

  .. code-block:: bash  

    sed -i 's/MODE=prod/MODE=library/g' Makefile 
    sed -i 's/COMP=gnu/COMP=intel/g' Makefile # - if using intel compiler

3. Link fftw to picsar (you may have to modify the following according to your machine)

  .. code-block:: bash  

    sed -i  's/FFTW3_LIB=\/usr\/lib\/x86_64-linux-gnu/FFTW3_LIB=$(FFTW_LIB_DIR)/g' Makefile
    sed -i  's/FFTW3_INCLUDE=\/usr\/include/FFTW3_INCLUDE=$(FFTW_INC_DIR)/g' Makefile

4. Install picsar as a library

  .. code-block:: bash  

    make lib
    export LIBPXR=$PWD/lib
    export LD_LIBRARY_PATH=${LIBPXR}:${LD_LIBRARY_PATH}
