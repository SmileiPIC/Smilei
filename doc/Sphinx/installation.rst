Install
-------

The requirements for installing and running :program:`Smilei` are:

* A C++ compiler (v4.8 or later recommended), optionally implementing openMP
* MPI libraries (*openmpi* recommended)
* HDF5 libraries compatible with your versions of C++ and MPI
* Python 2.7

Optional dependencies are:

* Doxygen
* Python modules: sphinx, h5py, numpy, matplotlib, pylab, pint
* ffmpeg

On a large cluster, refer to the administrator to install these requirements.
If you want to install :program:`Smilei` on your personal computer, refer to the following sections.

----

.. _installMac:

Install dependencies on Mac
^^^^^^^^^^^^^^^^^^^^^^^^^^^

Even if it possible to install all dependencies manually, we recommend a package manager software.


Via Macports
""""""""""""

This installation procedure relies on the software `MacPorts <https://www.macports.org>`_.

#. If you do not have it already, `install MacPorts <https://www.macports.org/install.php>`_.

   .. code-block:: bash

     sudo port -v selfupdate

#. In a terminal, run the following command to install the C++ compiler with MPI:
     
   .. code-block:: bash

     sudo port install openmpi-gcc48
     
   Then, to make this the default:
     
   .. code-block:: bash

     sudo port select --set mpi openmpi-gcc48-fortran
   
#. To install HDF5, run:
     
   .. code-block:: bash

     sudo port install hdf5 +gcc48+openmpi

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
   
   Then, to make this the default:
     
   .. code-block:: bash

     sudo port select --set python python27
     sudo port select --set python2 python27

#. in a new terminal window (to take into account of the above command) compile :program:`smilei` (see :ref:`compile`)

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

This installation procedure has been tested on OS X "El Capitan" 10.11.1

#. `HomeBrew <http://brew.sh>` does not need administrator privileges and can easily installed via:

   .. code-block:: bash

     ruby -e "$(curl -fsSL https://raw.githubusercontent.com/Homebrew/install/master/install)"

#. install the following packages using :program:`brew` to be able to compile and run :program:`smilei`

   .. code-block:: bash

     brew tap homebrew/science
     brew install gcc
     brew cask install java
     brew install makedepend
     HOMEBREW_CC=gcc-5 HOMEBREW_CXX=g++-5 brew install open-mpi --build-from-source
     brew install hdf5 --with-mpi
     brew install python

#. now you can compile :program:`smilei` (see :ref:`compile`)

#. install the following extra packages (in order of importance)

   .. code-block:: bash

     export LC_ALL=en_US.UTF-8
     export LANG=en_US.UTF-8
     pip install ipython h5py pint sphinx matplotlib pylab
     brew install doxygen



----

Install dependencies on Ubuntu
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
    
    Install the following packages from terminal

   .. code-block:: bash
     
     sudo apt-get install git openmpi-bin libhdf5-openmpi-dev build-essential python-dev




----

.. _compile:

Download and compile
^^^^^^^^^^^^^^^^^^^^

#. Download the latest tarball :ref:`here <latestVersion>`.

#. Extract the tarball at the location of your choice.
   Let us assume it is located in your home directory ``~/smilei/``.

#. In a terminal, go to that location and compile:
     
     .. code-block:: bash
       
       cd ~/smilei
       make

   To speedup un multiple CPUs:
     
     .. code-block:: bash
       
       make -j 4  # compile with 4 processors 
   
   Help on make alternatives:
     
     .. code-block:: bash
       
       make help

   examples:
     
     .. code-block:: bash
       
       make config=debug            # to have debugging output (slow)
       make config=noopenmp         # to deactivate OpenMP support
       make config="debug noopenmp" # to activate debugging without OpenMP
       make doc                     # to compile the documentation
   



