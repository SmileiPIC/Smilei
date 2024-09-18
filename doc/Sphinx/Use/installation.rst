Install
-------

Installing Smilei requires several steps:

#. Install compilers and libraries that Smilei needs (*dependencies*)
#. Download Smilei
#. Setup your environment (*environment variables*)
#. Compile

----

Install the dependencies
^^^^^^^^^^^^^^^^^^^^^^^^

The **necessary** dependencies are:

* A C++11 compiler, optionally implementing openMP version > 4.5.
* An MPI library (by default a version supporting ``MPI_THREAD_MULTIPLE``).
  IntelMPI or OpenMPI are recommended.
* The **parallel** HDF5 library compiled with your versions of C++ and MPI.
* Python 3+ with header files.

When compiling on GPU:

* The C++ compiler must be GPU-aware (typically ``nvc++`` for NVIDIA or ``clang`` for AMD)
* A CUDA or HIP compiler is necessary (typically ``nvcc`` for NVIDIA or ``hipcc`` for AMD)

Optional dependencies are:

* `Git <https://git-scm.com/>`_ for version control
* Python modules for post-processing: sphinx, h5py, numpy, matplotlib, pint
* `FFmpeg <https://ffmpeg.org/>`_ for converting animations to videos

There are various ways to install all dependencies, depending on the platform:

* :doc:`On MacOs<install_macos>`
* :doc:`On Linux<install_linux>`
* :doc:`On a supercomputer<install_supercomputer>`

If you have successfully installed these dependencies on other platforms,
please :doc:`contact us </Overview/partners>` and share!


----

Download the Smilei source
^^^^^^^^^^^^^^^^^^^^^^^^^^^

Clone the latest :program:`Smilei` version from Github:

.. code-block:: bash

  cd /path/of/your/choice/
  git clone https://github.com/SmileiPIC/Smilei.git

If you prefer a direct download, see :ref:`here <latestVersion>`.

----

Setup environment variables for compilation
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Several environment variables may be required, depending on your setup.

* ``SMILEICXX``: the MPI-C++ compiler.
  Defaults to ``mpicxx``.
* ``HDF5_ROOT_DIR``: the folder of the HDF5 library.
  Defaults to ``$HDF5_ROOT``.
* ``BUILD_DIR``: the folder where the compilation should occur.
  Defaults to ``./build``.
* ``PYTHONEXE``: the python executable to use in smilei.
  Defaults to ``python``.
* ``CXXFLAGS``: flags for the C++ compiler.
* ``LDFLAGS``: flags for the linker.
* ``GPU_COMPILER``: the compiler for CUDA or HIP (typically ``nvcc`` or ``hipcc``).
  Defaults to ``$CC``.
* ``GPU_COMPILER_FLAGS``: flags for ``$GPU_COMPILER``.

The command ``make help`` can give you some information about your environment.

----

.. _compile:

Compile Smilei
^^^^^^^^^^^^^^^^^^^^^^^^^^^

In a terminal, go to the folder where you downloaded :program:`Smilei` and use the commmand

.. code-block:: bash

  make

If the compilation is successful, you should now have a new ``smilei`` executable.

----

Advanced compilation options
^^^^^^^^^^^^^^^^^^^^^^^^^^^^

.. rubric:: Compile with several processors (fast compilation)

.. code-block:: bash

  make -j 4  # Compiles on 4 threads

.. rubric:: Compilation configuration with keyword "config"

.. code-block:: bash

  make config=noopenmp        # Without OpenMP support
  make config=no_mpi_tm       # Without a MPI library which supports MPI_THREAD_MULTIPLE
  make config=gpu_nvidia      # For Nvidia GPU acceleration
  make config=gpu_amd         # For AMD GPU acceleration
  make config=debug           # With debugging output (slow execution)
  make config=scalasca        # For the Scalasca profiler
  make config=advisor         # For Intel Advisor
  make config=vtune           # For Intel Vtune
  make config=inspector       # For Intel Inspector
  make config=detailed_timers # More detailed timers, but somewhat slower execution

It is possible to combine arguments above within quotes, for instance:

.. code-block:: bash

  make config="debug noopenmp" # With debugging output, without OpenMP

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

----

Compilation for GPU accelerated nodes
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

On GPU, two compilers are used: a C++ compiler for the main code
(defined by the variable ``$SMILEICXX``) and a compiler for 
``.cu`` CUDA files (defined by the variable ``$GPU_COMPILER``).
For NVIDIA, it is recommended to use the ``nvhpc`` software kit
which includes the compilers ``nvc++`` and ``nvcc``.
For AMD, the equivalent ``ROCm`` software kit includes ``clang`` and ``hipcc``.

Generally, several flags must be supplied to these compilers in order
to target properly your system architecture. They must
be supplied in ``$CXXFLAGS`` and ``$GPU_COMPILER_FLAGS``.
Please refer to the system administrators to find available compilers
and the required flags for your machine, as well as the commands
needed to load the correct environment.

The compilation of Smilei must include a special ``config`` keyword equal to either
``gpu_nvidia`` or ``gpu_amd``.
Two examples are provided as guidance:

.. code-block:: bash

  make -j 12 machine="jean_zay_gpu_A100" config="gpu_nvidia" # example for Nvidia GPU
  make -j 12 machine="adastra" config="gpu_amd"              # example for AMD GPU

In these cases, the environment variables were included in *machine files* that
you can find in ``scripts/compile_tools/machine/``.
Typically ``CXXFLAGS += -ta=tesla:cc80`` for ``nvhpc`` <23.4 and
``CXXFLAGS += -gpu=cc80 -acc`` for the more recent versions of ``nvhpc``.

.. warning::
  
  * The hdf5 module should be compiled with the nvidia/cray compiler;
    openmpi as well, but depending on the nvhpc module it might not be needed as it can be included in the nvhpc module.

For testing purposes, Smilei can be run on a linux PC with a good
GPU. As a guidance, read these :doc:`instructions for GPU on linux<install_linux_GPU>`.

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
