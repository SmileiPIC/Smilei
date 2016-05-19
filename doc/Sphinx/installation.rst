Install
-------

Before installing :program:`Smilei`, you need to install a few dependencies:

* A C++ compiler, optionally implementing openMP
* MPI libraries (*openmpi* recommended), supporting `MPI_THREAD_MULTIPLE`
* HDF5 libraries compatible with your versions of C++ and MPI
* Python 2.7

Optional dependencies are:

* Doxygen
* Python modules: sphinx, h5py, numpy, matplotlib, pylab, pint
* ffmpeg

On a large cluster, refer to the administrator to install these requirements.
If you want to install :program:`Smilei` on your personal computer, refer to the following sections.

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

This installation procedure has been tested on OS X "El Capitan" 10.11.1

#. `HomeBrew <http://brew.sh>`_ does not need administrator privileges and can easily installed via:

   .. code-block:: bash

     ruby -e "$(curl -fsSL https://raw.githubusercontent.com/Homebrew/install/master/install)"

#. install the following packages using :program:`brew` to be able to compile and run :program:`smilei`

   .. code-block:: bash

     brew tap homebrew/science
     brew cask install java
     brew install makedepend
     brew install  gcc5
     HOMEBREW_CC=gcc-5 HOMEBREW_CXX=g++-5 brew install -s openmpi --without-fortran --with-mpi-thread-multiple
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
    
Install the following packages from terminal:

  .. code-block:: bash
  
    sudo apt-get install git openmpi-bin libhdf5-openmpi-dev build-essential python-dev

:red:`Need details`


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

   .. rubric:: Compilation alternatives:
     
   .. code-block:: bash
     
     make -j 4                    # compile with 4 processors (fast)  
     make config=debug            # to have debugging output (slow)
     make config=noopenmp         # to deactivate OpenMP support
     make config="debug noopenmp" # to activate debugging without OpenMP
     make doc                     # to compile the documentation
     make help                    # to get some help on compilation
 
#. The next step is to :doc:`write a namelist <namelist>`.


