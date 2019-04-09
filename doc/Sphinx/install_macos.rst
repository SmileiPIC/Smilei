
Installing dependencies on Mac
---------------------------------

Even if it possible to install all dependencies manually, we recommend using a
package manager software.


Via Macports
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

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
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

This installation procedure has been tested on macOS 10.14.4

#. `HomeBrew <http://brew.sh>`_ can easily installed via:

   .. code-block:: bash

     /usr/bin/ruby -e "$(curl -fsSL https://raw.githubusercontent.com/Homebrew/install/master/install)"

#. Once installed, you need these packages:

   .. code-block:: bash

     brew install ipython iltommi/homebrew-brews/hdf5
     pip3 install h5py pint sphinx matplotlib scipy

#. Edit your ``.bash_profile`` hidden file located in your home folder:
   
   .. code-block:: bash

     open ~/.bash_profile
   
   and add the following lines at the end:
     
   .. code-block:: bash
  
     export OMPI_CXX=g++-8 
     export HDF5_ROOT_DIR=/usr/local
     export PYTHONEXE=python3

#. now you can compile :program:`smilei` (see :ref:`compile` for other options)

