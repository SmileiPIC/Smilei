
Prerequisites
---------------------------------

First, you will need to install Xcode and the Command Line Tools in order to be able to compile Smilei

   .. code-block:: bash

     xcode-select --install

and follow the instructions.


Brew : install Smilei
---------------------------------

This installation procedure has been tested on macOS 10.14.4

#. Install `HomeBrew <http://brew.sh>`_ via:

   .. code-block:: bash

     /usr/bin/ruby -e "$(curl -fsSL https://raw.githubusercontent.com/Homebrew/install/master/install)"

#. Install Smilei

   .. code-block:: bash

     brew install --HEAD iltommi/brews/smilei

   Smilei executables (``smilei`` and ``smilei_test``) and the python module are now accessible from everywhere.

#. Install python packages needed for the happi python module:

   .. code-block:: bash

     pip3 install ipython h5py pint sphinx matplotlib scipy
     

Documentation can be opened with

  .. code-block:: bash
  
     open /usr/local/opt/smilei/share/html/index.html     
     
To update Smilei with just type

  .. code-block:: bash
  
     brew upgrade --fetch-HEAD smilei

Brew : install dependencies
---------------------------------

In case you want ot keep a private version of Smilei where you can make changes to the core code, 
you might want to just install the Smilei dependencies to be able to compile Smilei from you directory:

#. install Smilei dependencies

   .. code-block:: bash
     
     brew install iltommi/brews/smilei --HEAD --only-dependencies

#. Edit your ``.bash_profile`` hidden file located in your home folder:
   
   .. code-block:: bash

     open ~/.bash_profile
   
   and add the following lines at the end:
     
   .. code-block:: bash
  
     export OMPI_CXX=g++-8 
     export HDF5_ROOT_DIR=/usr/local/opt/hdf5-parallel
     export PYTHONEXE=python3

#. In a new terminal window, you can now compile :program:`smilei` (see :ref:`compile` for other options)


Macports : install dependencies
---------------------------------

**Please note that these guidelines might be slightly outdated.**

If you find any error, please fill an issue on GitHub: https://github.com/SmileiPIC/Smilei/issues

This installation procedure relies on the software `MacPorts <https://www.macports.org>`_
that you can install following `these instructions <https://www.macports.org/install.php>`_.

#. In a terminal, run the following command to install the C++ compiler with MPI and HDF5:
     
   .. code-block:: bash

     sudo port install openmpi-gcc7 +threads
     sudo port select --set mpi openmpi-gcc7-fortran
     sudo port install hdf5 +openmpi+gcc7+threads
     
#. Edit your ``.bash_profile`` hidden file located in your home folder:
   
   .. code-block:: bash

     open ~/.bash_profile
   
   and add the following lines at the end:
     
   .. code-block:: bash

     export SMILEICXX=mpicxx
     export HDF5_ROOT_DIR=/opt/local/hdf5/lib/

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
