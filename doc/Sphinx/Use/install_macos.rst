
Install dependencies on MacOS
-----------------------------

#. Install `HomeBrew <http://brew.sh>`_ 

#. Install Smilei dependencies:

   .. code-block:: bash

     brew install python numpy hdf5-mpi

#. Install python packages needed for the ``happi`` python module:

   .. code-block:: bash

     pip3 install ipython h5py pint sphinx matplotlib scipy
     
#. Add the following lines to ``~/.bash_profile`` (or ``~/.zprofile`` of you're using zsh) :
     
   .. code-block:: bash
  
      export OMPI_CXX=g++-14
      export HDF5_ROOT_DIR=$(brew --prefix hdf5-mpi)

#. In a new terminal window, you can now compile :program:`smilei` (see :ref:`compile` for other options)



