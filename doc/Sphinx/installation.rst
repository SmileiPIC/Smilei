Install
-------

The requirements for installing and running :program:`Smilei` are:

* A C++ compiler (:red:`minimum version?`) with MPI (:red:`minimum version?`).
* The HDF5 libraries (:red:`minimum version?`) compatible with your version of MPI.
* Python 2.7.

Optional dependencies are:

* OpenMP (:red:`minimum version?`).
* The *Sphinx* documentation generator .
* Other python packages: h5py, numpy, matplotlib, pylab.
* ffmpeg.

On a large cluster, refer to the administrator to install these requirements.
If you want to install :program:`Smilei` on your personal computer, refer to the following sections.

----

.. _installMac:

Install dependencies on Mac
^^^^^^^^^^^^^^^^^^^^^^^^^^^

This installation procedure relies on the software `MacPorts <https://www.macports.org/>`_.
It is also possible to install all dependencies manually.

#. If you do not have it already, `install MacPorts <https://www.macports.org/install.php>`_.

#. In a terminal, run the following command to install the C++ compiler with MPI::
     
     $ sudo port install openmpi-gcc48
     
   Then, to make this the default::
     
     $ sudo port select --set mpi openmpi-gcc48-fortran
   
#. To install HDF5, run::
     
     $ sudo port install hdf5 +gcc48+openmpi
       
#. Optionally, to install openMP, 
   
   :red:`to do`
   
#. Edit your ``.bash_profile`` hidden file located in your home folder::
   
     $ open ~/.bash_profile
   
   and add the following lines in that file:
     
   .. code-block:: bash

     export SMILEICXX=mpicxx
     export HDF5_ROOT_DIR=/opt/local
     
   Depending on your system, you might need to use ``mpic++`` instead of ``mpicxx``.
  
#. Python should be already installed by default, but in case you need
   a specific version, run::
   
     $ sudo port install python27
   
   Then, to make this the default::
     
     $ sudo port select --set python python27
     $ sudo port select --set python2 python27

#. If you wish to run the Python post-processing scripts provided in :program:`Smilei`,
   you need several modules (numpy, matplotlib, pylab, h5py). We recommend to install
   :program:`IPython` which includes some of these::
   
     $ sudo port install py27-ipython
   
   Then, for h5py::
     
     $ sudo port install py27-h5py
   
#. If you need to build the documentation as well, refer to the `README` provided
   in :program:`Smilei`.


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
   * ``make -j4`` to compile with 4 processors.
   * ``make doc`` to compile the documentation.



