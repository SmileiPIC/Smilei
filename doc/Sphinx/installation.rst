Requirements and installation
-----------------------------

The requirements for installing and running :program:`Smilei` are:

* A C++ compiler with MPI.
* The HDF5 libraries compatible with your version of MPI.
* Python 2.7.

Optional dependencies are:

* OpenMP.
* The *Sphinx* documentation generator.
* Other python packages: h5py, numpy, matplotlib, pylab.
* ffmpeg.

The installation steps depend on your system.
On a large cluster, refer to the administrator to install the requirements above.
If you want to install :program:`Smilei` on your personal computer, refer to the section
:ref:`Mac <installMac>` or :ref:`Ubuntu <installUbuntu>`.


----

Downloading and compiling
^^^^^^^^^^^^^^^^^^^^^^^^^

Before you compile :program:`Smilei`, make sure you have the dependencies stated above.

#. Download the latest tarball :ref:`here <latestVersion>`.

#. Extract the tarball at the location of your choice.
   Here we suppose it is located in your home directory ``~/Smilei/``.

#. :red:`setting PATH or other environment variables ?`

#. In a terminal, go to that location and compile::
     
     $ cd ~/Smilei
     $ make
   
   Alternates:
     
   * ``make debug`` to have debugging output (slow).
   * ``make -j4`` to compile with 4 processors.
   * ``make doc`` to compile the documentation.


----

Running
^^^^^^^

.. rst-class:: inprogress
  
  In progress ...



----

.. _installMac:

Install dependencies on Mac
^^^^^^^^^^^^^^^^^^^^^^^^^^^

This installation procedure relies on the software `MacPorts <https://www.macports.org/>`_.
It is also possible to install all dependencies manually.

#. If you do not have it already, `install MacPorts <https://www.macports.org/install.php>`_.

#. In a terminal, run the following command to install the C++ compiler with MPI::
     
     $ sudo port install openmpi-gcc48
   
   :red:`to be confirmed` ... 
   
   Do **not** use ``mpich`` instead of ``openmpi``.
   
   :red:`mpi-select or something like that ?` ... 

   It may take a while.
   
   You may change the version of GCC but it is not guaranteed to work with MPI or HDF5.

#. To install HDF5, run::
     
     $ sudo port install hdf5 +gcc48+openmpi
  
   :red:`to be confirmed` ... 
   
   :red:`select a variant ?` ... 
   
#. Edit your ``.bash_profile`` hidden file located in your home folder::
   
     $ open ~/.bash_profile
   
   and add the following lines in that file:
     
   .. code-block:: bash

     export SMILEICXX=mpicxx # This might require another executable
     export HDF5_ROOT_DIR=/opt/local     
  
#. Python should be already installed by default, but in case you need
   a specific version, run::
   
     $ sudo port install python27
   
   (follow the instructions on screen to make this version default).
   
#. If you wish to run the Python post-processing scripts provided in :program:`Smilei`,
   you need several modules (numpy, matplotlib, pylab, h5py). We recommend to install
   :program:`IPython` which includes some of these::
   
     $ sudo port install py27-ipython
   
   Then, for h5py::
     
     $ sudo port install py27-h5py
   
#. If you need to build the documentation as well, refer to the `README` provided
   in :program:`Smilei`.


----

.. _installUbuntu:

Install dependencies on Ubuntu
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^




