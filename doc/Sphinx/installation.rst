Install
-------

Before installing :program:`Smilei`, you need to install a few dependencies:

* A C++11 compiler, optionally implementing openMP
* MPI libraries supporting ``MPI_THREAD_MULTIPLE``
* HDF5 libraries compatible with your versions of C++ and MPI
* Python 2.7 or Python 3 (with header files)

Optional dependencies are:

* Git
* Python modules: sphinx, h5py, numpy, matplotlib, pylab, pint
* ffmpeg
* the `Picsar <http://picsar.net>`_ library

On a large cluster, refer to the administrator to install these requirements.
If you want to install :program:`Smilei` on your personal computer, refer to the following sections.

To get some help on compilation and the environment variables you can change in order
to have a successful compilation, you can type ``make help``.

----

Install dependencies on Mac
^^^^^^^^^^^^^^^^^^^^^^^^^^^

Even if it possible to install all dependencies manually, we recommend using a
package manager software.


Via Macports
""""""""""""

This installation procedure relies on the software `MacPorts <https://www.macports.org>`_
that you can install following `these instructions <https://www.macports.org/install.php>`_.

#. In a terminal, run the following command to install the C++ compiler with MPI and HDF5:

   .. code-block:: bash

     sudo port install openmpi-gcc7 +threads
     sudo port select --set mpi openmpi-gcc7-fortran
     sudo port install hdf5-18 +openmpi+gcc7+threads

#. Edit your ``.bash_profile`` hidden file located in your home folder:

   .. code-block:: bash

     open ~/.bash_profile

   and add the following lines at the end:

   .. code-block:: bash

     export SMILEICXX=mpicxx
     export HDF5_ROOT_DIR=/opt/local/hdf5-18/lib/

#. Python should be already installed by default, but in case you need
   a specific version, run:

   .. code-block:: bash

     sudo port install python27
     sudo port select --set python python27
     sudo port select --set python2 python27

#. If you wish to run the Python post-processing scripts provided in :program:`Smilei`,
   you need several modules (numpy, matplotlib, pylab, h5py, sphinx, pint).
   We recommend to install :program:`IPython` which includes some of these.

   .. code-block:: bash

     sudo port install py27-ipython # nicer python console
     sudo port install py27-h5py    # mandatory for opening any HDF5 file
     sudo port install py27-pint    # only for auto unit conversion
     sudo port install py27-sphinx  # only for building the doc


Via HomeBrew
""""""""""""

This installation procedure has been tested on macOS 10.13

#. `HomeBrew <http://brew.sh>`_ can easily installed via:

   .. code-block:: bash

     ruby -e "$(curl -fsSL https://raw.githubusercontent.com/Homebrew/install/master/install)"

#. Once installed, you need these packages:

   .. code-block:: bash

     brew install gcc python numpy
     brew install openmpi --with-mpi-thread-multiple
     brew install hdf5 --with-mpi
     export LC_ALL=en_US.UTF-8
     export LANG=en_US.UTF-8
     pip3 install ipython h5py pint sphinx matplotlib scipy pylab

#. To be able to use the gcc with openmpi, you need to set the ``OMPI_CXX`` variable :

   .. code-block:: bash

     export OMPI_CXX=g++-7 # the number version might vary

#. You can put the above line in a shell rc file (e.g. ``.bash_profile``)
   or you can just add it before the ``make`` command (``OMPI_CXX=g++-7 make`` ...)

#. now you can compile :program:`smilei` (see :ref:`compile` for other options)


----

Install dependencies on Linux
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Fedora
""""""

.. code-block:: bash

  dnf install -y gcc-c++ hdf5-openmpi hdf5-openmpi-devel openmpi-devel git which findutils python python-devel
  dnf install -y h5py ipython python2-pint sphinx python2-matplotlib



Debian (Ubuntu, Mint etc...)
""""""""""""""""""""""""""""

1. Install these packages

  .. code-block:: bash

    sudo apt-get install python-h5py ipython python-pint python-sphinx python-matplotlib python2-dev  pyhton-numpy

Since the system ``openmpi`` is not compiled with ``--enable-mpi-thread-multiple``, a manual installation is required :

2. Choose a path whet to install dependencies by setting the environment variable ``INSTALL_DIR``. e.g. :

  .. code-block:: bash

    export INSTALL_DIR=/usr/local


3. Download `OpenMPI <https://www.open-mpi.org/software/ompi>`_ and install.
   You may choose any ``${INSTALL_DIR}``.

  .. code-block:: bash

    tar zxvf openmpi-1.10.2.tar.gz
    cd openmpi-1.10.2
    ./configure --prefix=${INSTALL_DIR}/openmpi --enable-mpi-thread-multiple --enable-mpirun-prefix-by-default
    make
    make install

  Set environment variables in your `~/.bashrc` or `~/.bash_profile` file.

  .. code-block:: bash

    export PATH=${INSTALL_DIR}/openmpi/bin:${PATH}
    export LD_LIBRARY_PATH=${INSTALL_DIR}/openmpi/lib:${LD_LIBRARY_PATH}


4. Download `HDF5 <https://support.hdfgroup.org/HDF5>`_ and install

  .. code-block:: bash

    tar zxvf hdf5-1.8.16.tar.gz
    cd hdf5-1.8.16
    ./configure --prefix=${INSTALL_DIR}/hdf5 --enable-parallel --with-pic --enable-linux-lfs --enable-shared --enable-production=yes --disable-sharedlib-rpath --enable-static CC=mpicc FC=mpif90
    make
    make install

  Set environment variables in your `~/.bashrc` or `~/.bash_profile` file.

  .. code-block:: bash

    export PATH=${INSTALL_DIR}/hdf5/bin:${PATH}
    export LD_LIBRARY_PATH=${INSTALL_DIR}/hdf5/lib:${LD_LIBRARY_PATH}
    export HDF5_ROOT_DIR=${INSTALL_DIR}/hdf5


----


The optional Picsar library
^^^^^^^^^^^^^^^^^^^^^^^^^^^

A Pseudo-Spectral Analytical Time Domain solver (PSATD) for Maxwell equations
is available for use in Smilei via the `Picsar <http://picsar.net>`_ library.
The PSATD solver is an FFT-based high-order Maxwell solver.
Picsar uses the `FFTW <http://www.fftw.org>`_ library.

Install FFTW
""""""""""""

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
"""""""""""""""""""""""""""

1. Download the `latest version of Picsar <git@bitbucket.org:berkeleylab/picsar.git>`_

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


----

Install dependencies on other systems
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

If you have successfully installed these dependencies on other platforms, please
:doc:`contact us <partners>` and share!

----

.. _compile:

Download and compile Smilei
^^^^^^^^^^^^^^^^^^^^^^^^^^^

#. Download the latest :program:`Smilei` tarball :ref:`here <latestVersion>`.

#. Extract the tarball at the location of your choice.
   Let us assume it is located in your home directory ``~/smilei/``.

#. In a terminal, go to that location and compile:

   .. code-block:: bash

     cd ~/smilei
     make

#. The next step is to :doc:`write a namelist <namelist>`.

----

Advanced compilation options
^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Several ``make`` options are available:

.. code-block:: bash

  make -j 4                    # Compiles with 4 procs (fast compilation)
  make config=debug            # With debugging output (slow execution)
  make config=noopenmp         # Without OpenMP support
  make config="debug noopenmp" # With debugging output, without OpenMP
  make print-XXX               # Prints the value of makefile variable XXX
  make env                     # Prints the values of all makefile variables
  make help                    # Gets some help on compilation
  sed -i 's/PICSAR=FALSE/PICSAR=TRUE/g' makefile; make -j4 #To enable calls for PSATD solver from picsar


Each machine may require a specific configuration (environment variables, modules, etc.).
Such instructions may be included, from a file of your choice, via the ``machine`` argument:

.. code-block:: bash

  make machine=my_machine_file

where ``my_machine_file`` is a file, located in ``scripts/CompileTools/machine``, containing
the lines of command to be executed before compilation.

If you successfully write such a file for a common supercomputer, please share it
with developpers so that it can be included in the next release of :program:`Smilei`.

----

Advanced compilation options for profiling/tracing tools
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Several ``make`` options are available in order to compile Smilei ready for
specific performance analysis and debugging tools.

Smilei performance can be analyzed with :
- Scalasca
- Intel Advisor
- Intel Vtune
- Intel Inspector

.. code-block:: bash

  make config="scalasca"             : compilation for scalasca (required the Scalasca profiler)
  make config="advisor"              : compilation for Intel Advisor (required the Intel suite)
  make config="vtune"                : compilation for Intel Vtune (required the Intel suite)
  make config="inspector"            : compilation for Intel Inspector (required the Intel suite)
----

Advanced compilation options for detailed timers
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

The code contains more timers than in the default configuration that can
be activated using a specific option called ``detailed_timers``.
This flag will add the following line to the compilation flags ``-D__DETAILED_TIMERS``
and therefore activate in the code their computation.
These timers are not available in the default compilation
because they are more intrusive and may impact the overall performance in production.

Some of the timers are called inside patches for specific operators
such as the particle pusher. The code first computes the average over the patches
for all MPI domain before computing the final min, max and mean values between
MPI processes.

The final information is therefore an average that does not reflect the
load imbalance between patches of a MPI domain.

the following ``make`` command will activate the additional timers:

.. code-block:: bash

  make config="detailed_timers"      : compilation that activates the detailed timers inside the patch

Additional timers will be shown at the end of the simulation and are also
in ``profile.txt``

----

Compile the documentation
^^^^^^^^^^^^^^^^^^^^^^^^^

If you have installed the python module ``sphinx``, you can compile the documentation
(which you are currently reading) with:

.. code-block:: bash

   make doc

This creates a local *html* website accessible in your ``build/html/`` folder.

----

.. _installModule:

Install the happi module
^^^^^^^^^^^^^^^^^^^^^^^^

A python module, ``happi``, is provided to view, extract and post-process data from
all the diagnostics.
There are several ways to load this module in python.

1. Recommended:

  .. code-block:: bash

    make happi

  This has to be done only once, unless you move the smilei directory elsewhere.
  This command creates a small file in the Python *user-site* directory that tells python
  where to find the module.
  To remove it use the command ``make uninstall_happi``.

  The module will directly be accessible from *python*::

    >>> import happi

2. Alternative: Execute the ``Diagnostics.py`` script from python

  Adding a new *python* module is not always possible.
  Instead, we provide the script ``Diagnostics.py`` which is able to find the ``happi``
  module and import it into *python*.

  You may add the following command in your own python script::

    >>> execfile("/path/to/Smilei/scripts/Diagnostics.py")
