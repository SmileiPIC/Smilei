Install
-------

The requirements for installing and running :program:`Smilei` are:

* A C++ compiler (v4.8 or later recommended) with MPI (*openmpi* is preferred)
* The HDF5 libraries compatible with your versions of C++ and MPI
* Python 2.7

Optional dependencies are:

* OpenMP
* Doxygen
* Python modules: sphinx, h5py, numpy, matplotlib, pylab, pint
* ffmpeg

On a large cluster, refer to the administrator to install these requirements.
If you want to install :program:`Smilei` on your personal computer, refer to the following sections.

----

.. _installMac:

Install dependencies on Mac
^^^^^^^^^^^^^^^^^^^^^^^^^^^

This installation procedure relies on the software `MacPorts <https://www.macports.org/>`_.
It is also possible to install all dependencies manually.

#. If you do not have it already, `install MacPorts <https://www.macports.org/install.php>`_.

#. In a terminal, run the following command to install the C++ compiler with MPI:
     
   .. code-block:: bash

     $ sudo port install openmpi-gcc48
     
   Then, to make this the default:
     
   .. code-block:: bash

     $ sudo port select --set mpi openmpi-gcc48-fortran
   
#. To install HDF5, run:
     
   .. code-block:: bash

     $ sudo port install hdf5 +gcc48+openmpi
       
#. Optionally, to install openMP
   
   :red:`to do`
   
#. Edit your ``.bash_profile`` hidden file located in your home folder:
   
   .. code-block:: bash

     $ open ~/.bash_profile
   
   and add the following lines in that file:
     
   .. code-block:: bash

     export SMILEICXX=mpicxx
     export HDF5_ROOT_DIR=/opt/local
     
   Depending on your system, you might need to use ``mpic++`` instead of ``mpicxx``.
   Default for `SMILEICXX`  is `mpic++`, `HDF5_ROOT_DIR` is empty.

#. Python should be already installed by default, but in case you need
   a specific version, run:
   
   .. code-block:: bash

     $ sudo port install python27
   
   Then, to make this the default:
     
   .. code-block:: bash

     $ sudo port select --set python python27
     $ sudo port select --set python2 python27

#. If you wish to run the Python post-processing scripts provided in :program:`Smilei`,
   you need several modules (numpy, matplotlib, pylab, h5py, sphinx, pint).
   We recommend to install :program:`IPython` which includes some of these.
   
   .. code-block:: bash

     $ sudo port install py27-ipython
     $ sudo port install py27-h5py
     $ sudo port install py27-sphinx # only for building the doc
     $ sudo port install py27-pint   # only for auto unit conversion


----

Install dependencies on Ubuntu
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

.. rst-class:: inprogress
  
  In progress ...




----

Download and compile
^^^^^^^^^^^^^^^^^^^^

#. Download the latest tarball :ref:`here <latestVersion>`.

#. Extract the tarball at the location of your choice.
   Let us assume it is located in your home directory ``~/Smilei/``.

#. In a terminal, go to that location and compile::
     
     $ cd ~/Smilei
     $ make
   
   Alternates:
     
   * ``make debug`` to have debugging output (slow).
   * ``make openmp`` to activate OpenMP support
   * ``make -j4`` to compile with 4 processors.
   * ``make doc`` to compile the documentation.
   



