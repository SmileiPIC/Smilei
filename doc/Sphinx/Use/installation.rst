Install
-------

Before installing :program:`Smilei`, you need to install a few dependencies:

* A C++11 compiler, optionally implementing openMP version > 4.5
  (gcc users: v6.0 or newer recommended)
* an MPI library (by default a version supporting ``MPI_THREAD_MULTIPLE``
  is required: v4.0 or newer recommended)
* an HDF5 library compatible with your versions of C++ and MPI
* Python 2.7 or Python 3+ (with header files)

Optional dependencies are:

* Git
* Python modules: sphinx, h5py, numpy, matplotlib, pint
* ffmpeg
* CUDA for NVIDIA GPUs or HIP-SYCL for AMD GPUs (it is recommended to use the already installed software stack and the support team of a supercomputer you have access to). 

----

Install the dependencies
^^^^^^^^^^^^^^^^^^^^^^^^

There are various ways to install all dependencies, depending on the platform:

* :doc:`On MacOs<install_macos>`
* :doc:`On Linux<install_linux>`
* :doc:`On a supercomputer<install_supercomputer>`

The command ``make help`` can give you some information about your environment.

If you have successfully installed these dependencies on other platforms,
please :doc:`contact us </Overview/partners>` and share!


----

Setup environment variables for compilation
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Several environment variables may be required, depending on your setup.

* ``SMILEICXX``: the MPI-C++ compiler.
  Defaults to ``mpicxx``.
* ``HDF5_ROOT_DIR``: the folder for the HDF5 library.
  Defaults to ``$HDF5_ROOT``.
* ``BUILD_DIR``: the folder where the compilation should occur.
  Defaults to ``./build``.
* ``PYTHONEXE``: the python executable to use in smilei.
  Defaults to ``python``.

The usual ``CXXFLAGS`` and ``LDFLAGS`` can also be used to pass other
arguments to the compiler and linker.


----

.. _compile:

Download and compile
^^^^^^^^^^^^^^^^^^^^^^^^^^^

#. Clone the latest :program:`Smilei` version from Github:

   .. code-block:: bash
    
     cd /path/of/your/choice/
     git clone https://github.com/SmileiPIC/Smilei.git
    
   If you do not have ``git``, you can dowload a tarball :ref:`here <latestVersion>`
   and extract it in a new folder.

#. In a terminal, go to that location and compile:

   .. code-block:: bash

     cd Smilei
     make
   
   If the compilation is successful, you should now have a new ``smilei`` executable.

#. The next step is to :doc:`write a namelist <namelist>`.

----

Advanced compilation options
^^^^^^^^^^^^^^^^^^^^^^^^^^^^

.. rubric:: Compile with several processors (fast compilation)

.. code-block:: bash

  make -j 4

.. rubric:: Compilation configuration with keyword "config"

.. code-block:: bash

  make config=debug                        # With debugging output (slow execution)
  make config=noopenmp                     # Without OpenMP support
  make config=no_mpi_tm                    # Without a MPI library which supports MPI_THREAD_MULTIPLE
  make config=scalasca                     # For the Scalasca profiler
  make config=advisor                      # For Intel Advisor
  make config=vtune                        # For Intel Vtune
  make config=inspector                    # For Intel Inspector
  make config=detailed_timers              # More detailed timers, but somewhat slower execution
  make config=omptasks                     # use OpenMP task parallelization, not supported by old compilers
  make config=part_event_tracing_tasks_off # trace the use particle operators, without task parallelization
  make config=part_event_tracing_tasks_on  # trace the use particle operators, with OpenMP task parallelization
  make config="gpu_nvidia noopenmp"        # For Nvidia GPU acceleration
  make config="gpu_amd"                    # For AMD GPU acceleration

It is possible to combine arguments above within quotes, for instance:

.. code-block:: bash

  make config="debug noopenmp" # With debugging output, without OpenMP

However, some arguments may not be compatible, e.g. ``noopenmp`` and ``omptasks``. 

.. rubric:: Obtain some information about the compilation

.. code-block:: bash

  make print-XXX               # Prints the value of makefile variable XXX
  make env                     # Prints the values of all makefile variables
  make help                    # Gets some help on compilation

.. rubric:: Machine-specific compilation

Each machine may require a specific configuration (environment variables,
modules, etc.). These instructions may be included in a file of your choice,
via the ``machine`` argument:

.. code-block:: bash

  make machine=my_machine_file

where ``my_machine_file`` is a file, located in
``scripts/compile_tools/machine``, containing the lines of command to be
executed before compilation. If you successfully write such a file for
a common supercomputer, please share it with developpers so that it can
be included in the next release of :program:`Smilei`.


.. rubric:: Compilation for GPU accelerated nodes:

As each supercomputer has a different environnment to compile for GPUs and since the nvhpc + CUDA/ cray + HIP modules evolve quickly, a machine file is required for the compilation.
Several machine files are already available as an example in smilei/scripts/compile_tools/machine/ ; such as: jean_zay_gpu_V100, jean_zay_gpu_A100, adastra, ruche_gpu2.

Typically we need it to specify ACCELERATOR_GPU_FLAGS += -ta=tesla:cc80 for nvhpc <23.4 and ACCELERATOR_GPU_FLAGS += -gpu=cc80 -acc for the more recent versions of nvhpc.

.. code-block:: bash

	make -j 12 machine="jean_zay_gpu_A100" config="gpu_nvidia noopenmp verbose" # for Nvidia GPU
	make -j 12 machine="adastra" config="gpu_amd" 			            # for AMD GPU


Furthermore, here are 2 examples of known working ennvironments, first for AMD GPUs, second for Nvidia GPUs:

.. code-block:: bash

	module purge
	module load craype-accel-amd-gfx90a craype-x86-trento
	module load PrgEnv-cray/8.3.3
	module load cpe/23.02
	module load cray-mpich/8.1.24 cray-hdf5-parallel/1.12.2.1 cray-python/3.9.13.1
	module load amd-mixed/5.2.3

.. code-block:: bash

	module purge
	module load anaconda-py3/2020.11  # python is fine as well if you can pip install the required modules
	module load nvidia-compilers/23.1
	module load cuda/11.2
	module load openmpi/4.1.1-cuda
	module load hdf5/1.12.0-mpi-cuda
	# For HDF5, note that module show can give you the right path
	export HDF5_ROOT_DIR=/DIRECTORY_NAME/hdf5/1.12.0/pgi-20.4-HASH/

Note: 

* we are aware of issues with CUDA >12.0, fixes are being tested but are not deployed yet. We recommend CUDA 11.x at the moment.
* The hdf5 module should be compiled with the nvidia/cray compiler ; openmpi as well, but depending on the nvhpc module it might not be needed as it can be included in the nvhpc module 

----

.. _vectorization_flags:

Optimization and vectorization options explained
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

To tune optimization and vectorization options, :program:`Smilei` uses the *machine files* described above. They contain compiler options for specific hardware architectures or processor families.

This :doc:`page <optimization_flags>` explains in detail optimization flags used in machine files and therefore how to generate your own machine file.

----

Create the documentation
^^^^^^^^^^^^^^^^^^^^^^^^

If you have installed the python module ``sphinx``, you can create the
documentation (which you are currently reading) with:

.. code-block:: bash

   make doc

This creates a local *html* website accessible in your ``build/html/`` folder.

----

.. _installModule:

Install the happi module
^^^^^^^^^^^^^^^^^^^^^^^^

A python module, ``happi``, is provided to view, extract and post-process
data from all the diagnostics.
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
  Instead, we provide the script ``Diagnostics.py`` which is able to find the
  ``happi`` module and import it into *python*.

  You may add the following command in your own python script::

    >>> execfile("/path/to/Smilei/scripts/Diagnostics.py")

----

Install the ``smilei_tables`` tool
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Generation of the tables is handled by an external tools.
A full documentation is available on :doc:`the dedicated page <tables>`.
