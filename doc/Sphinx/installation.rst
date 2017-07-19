Install
-------

Before installing :program:`Smilei`, you need to install a few dependencies:

* A C++11 compiler, optionally implementing openMP
* MPI libraries (**openmpi** recommended), supporting ``MPI_THREAD_MULTIPLE``
* HDF5 libraries compatible with your versions of C++ and MPI
* Python 2.7 (with header files)
* ``make``

Optional dependencies are:

* Git
* Doxygen
* Python modules: sphinx, h5py, numpy, matplotlib, pylab, pint
* ffmpeg

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

     sudo port install openmpi-gcc5 +threads
     sudo port select --set mpi openmpi-gcc5-fortran
     sudo port install hdf5 +openmpi+gcc5+threads
     
#. Edit your ``.bash_profile`` hidden file located in your home folder:
   
   .. code-block:: bash

     open ~/.bash_profile
   
   and add the following lines at the end:
     
   .. code-block:: bash

     export SMILEICXX=mpicxx
     export HDF5_ROOT_DIR=/opt/local
     
   Depending on your system, you might need to use ``mpic++`` instead of ``mpicxx``.

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
     sudo port install doxygen      # only for building the reference C++ doc


Via HomeBrew
""""""""""""

This installation procedure has been tested on OS X 10.12

#. `HomeBrew <http://brew.sh>`_ does not need administrator privileges and can easily installed via:

   .. code-block:: bash

     ruby -e "$(curl -fsSL https://raw.githubusercontent.com/Homebrew/install/master/install)"

#. install the following packages using :program:`brew` to be able to compile and run :program:`smilei`

   .. code-block:: bash

     brew tap homebrew/science
     brew install gcc
     brew install openmpi --with-mpi-thread-multiple
     brew install hdf5 --with-mpi     
     brew install python numpy

#. Now you need to set the ``OMPI_CXX`` to the homebrew ``g++`` (``g++-7`` or similar):
     
   .. code-block:: bash

     export OMPI_CXX=g++-7

#. Alternatively you can put this line variable in a shell rc file (e.g. ``.bash_profile``) 
   or you can just add it before the ``make`` command (``OMPI_CXX=g++-7 make`` ...)

#. now you can compile :program:`smilei` (see :ref:`compile`)

#. install the following extra packages (in order of importance)

   .. code-block:: bash

     export LC_ALL=en_US.UTF-8
     export LANG=en_US.UTF-8
     pip install ipython h5py pint sphinx matplotlib pylab



----

Install dependencies on Ubuntu
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
    
On Ubuntu 16.04
"""""""""""""""

Install the following packages from terminal:

  .. code-block:: bash
  
    sudo apt-get install git openmpi-bin libhdf5-openmpi-dev build-essential python-dev

On older release
""""""""""""""""

A manual installation is required :

1. Download `OpenMPI <https://www.open-mpi.org/software/ompi>`_

  .. code-block:: bash
  
    $ taz zxvf openmpi-1.10.2.tar.gz
    $ cd openmpi-1.10.2
    $ ./configure --prefix=${INSTALL_DIR}/openmpi-1.10.2 --enable-mpi-thread-multiple --enable-mpirun-prefix-by-default
    $ make
    $ make install
    $ export PATH=${INSTALL_DIR}/openmpi-1.10.2/bin:${PATH}
    $ export LD_LIBRARY_PATH=${INSTALL_DIR}/openmpi-1.10.2/lib:${LD_LIBRARY_PATH}


2. Download `HDF5 <https://support.hdfgroup.org/HDF5>`_

  .. code-block:: bash
  
    $ tar zxvf hdf5-1.8.16.tar.gz
    $ cd hdf5-1.8.16
    $ ./configure --prefix=${INSTALL_DIR}/hdf5-1.8.16 --enable-parallel --with-pic --enable-linux-lfs --enable-shared --enable-production=yes --disable-sharedlib-rpath --enable-static CC=mpicc FC=mpif90
    $ make
    $ make install
    $ export PATH=${INSTALL_DIR}/hdf5-1.8.16/bin:${PATH}
    $ export LD_LIBRARY_PATH ${INSTALL_DIR}/hdf5-1.8.16/lib:${LD_LIBRARY_PATH}
    $ # set HDF5 variable used in SMILEI makefile
    $ export HDF5_ROOT_DIR=${INSTALL_DIR}/hdf5-1.8.16


----

Install dependencies on other systems
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

If you have successfully installed these dependencies on other platforms, please
:doc:`contact us <partners>` and share!

----

.. _compile:

Download and compile
^^^^^^^^^^^^^^^^^^^^

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


Each machine may require a specific configuration (environment variables, modules, etc.).
Such instructions may be included, from a file of your choice, via the ``machine`` argument:

.. code-block:: bash
  
  make machine=my_machine_file

where ``my_machine_file`` is a file, located in ``scripts/CompileTools/machine``, containing
the lines of command to be executed before compilation.

If you successfully write such a file for a common supercomputer, please share it
with developpers so that it can be included in the next release of :program:`Smilei`.
 


----

Compile the documentation
^^^^^^^^^^^^^^^^^^^^^^^^^

There are two types of documentation:

#. `Doxygen`, only useful for developers
#. `Sphinx`, for users: the one you are currently reading

They can be compiled with the following commands:

.. code-block:: bash

   make doc     # Compiles all the documentation
   make doxygen # Compiles only the `doxygen` doc
   make sphinx  # Compiles only the `sphinx` doc


----

.. _installModule:

Install the Smilei python module
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

A python module is provided to view, extract and post-process data from all the diagnostics.
There are several ways to load this module in python.

1. Recommended: Install Smilei's module
  
  .. code-block:: bash
    
    make install_python
  
  This has to be done only once, unless you move the smilei directory elsewhere.
  This command creates a small file in the Python `user-site` directory that tells python
  where to find the module.
  To remove it use the command ``make uninstall_python``.
  
  The module is called ``Smilei`` and will directly be accessible from *python*::
    
    from Smilei import *

2. Alternative: Execute the ``Diagnostics.py`` script from python 
  
  Adding a new *python* module is not always possible.
  Instead, we provide the script ``Diagnostics.py`` which is able to find the ``Smilei``
  module and import it into *python*.
  
  You may add the following command in your own python script::
  
    execfile("/path/to/Smilei/scripts/Diagnostics.py")



