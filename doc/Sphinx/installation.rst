Install
-------

Before installing :program:`Smilei`, you need to install a few dependencies:

* A C++11 compiler, optionally implementing openMP
* an MPI library (by default a version supporting ``MPI_THREAD_MULTIPLE`` is required)
* an HDF5 library compatible with your versions of C++ and MPI
* Python 2.7 or Python 3 (with header files)

Optional dependencies are:

* Git
* Python modules: sphinx, h5py, numpy, matplotlib, pylab, pint
* ffmpeg
* the `Picsar <http://picsar.net>`_ library: see :doc:`this documentation<install_PICSAR>`

----

Install the dependencies
^^^^^^^^^^^^^^^^^^^^^^^^

There are various ways to install all dependencies, depending on the platform:

* :doc:`On MacOs<install_macos>`
* :doc:`On Linux<install_linux>`
* :doc:`On a supercomputer<install_supercomputer>`

The command ``make help`` can give you some information about your environment.

If you have successfully installed these dependencies on other platforms, please
:doc:`contact us <partners>` and share!

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

     cd ~/Smilei
     make
   
   If the compilation is successful, you should now have a new ``smilei`` executable.

#. The next step is to :doc:`write a namelist <namelist>`.

----

Advanced compilation options
^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Several ``make`` options are available:

.. code-block:: bash

  make -j 4                    # Compiles with 4 procs (fast compilation)
  make config=debug            # With debugging output (slow execution)
  make config=noopenmp         # Without OpenMP support
  make config="debug noopenmp" # With debugging output, without OpenMP
  make config=no_mpi_tm        # Without a MPI library which supports MPI_THREAD_MULTIPLE
  make print-XXX               # Prints the value of makefile variable XXX
  make env                     # Prints the values of all makefile variables
  make help                    # Gets some help on compilation
  sed -i 's/PICSAR=FALSE/PICSAR=TRUE/g' makefile; make -j4 #To enable calls for PSATD solver from picsar


----

.. _vectorization_flags:

Options for SIMD vectorization
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

The :doc:`SIMD vectorization <vectorization>` of :program:`Smilei` uses ``#pragma omp simd``.
To be enabled, you must provide appropriate options to your compiler through
the environment variable ``CXXFLAGS``.

For instance, :program:`Smilei` has been tested on
Intel processors (Skylake 8168) with an Intel environment.
The following flags provide a good performance:

.. code-block:: bash
  
  -xCOMMON-AVX512 -ip -ipo -inline-factor=1000 -D__INTEL_SKYLAKE_8168

The vectorization must also be activated :ref:`in the namelist <Vectorization>`.

----

Machine-specific compilation
^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Each machine may require a specific configuration (environment variables, modules, etc.).
Such instructions may be included, from a file of your choice, via the ``machine`` argument:

.. code-block:: bash

  make machine=my_machine_file

where ``my_machine_file`` is a file, located in ``scripts/CompileTools/machine``, containing
the lines of command to be executed before compilation.

If you successfully write such a file for a common supercomputer, please share it
with developpers so that it can be included in the next release of :program:`Smilei`.

----

Compilation options for profiling/tracing tools
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Several ``make`` options are available in order to compile Smilei ready for
specific performance analysis and debugging tools:

- Scalasca
- Intel Advisor
- Intel Vtune
- Intel Inspector

.. code-block:: bash

  make config="scalasca"   # compilation for scalasca (required the Scalasca profiler)
  make config="advisor"    # compilation for Intel Advisor (required the Intel suite)
  make config="vtune"      # compilation for Intel Vtune (required the Intel suite)
  make config="inspector"  # compilation for Intel Inspector (required the Intel suite)

----

Compilation options for detailed timers
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

The code contains optional timers for a more detailed timing
characterization. They are more intrusive
and may impact the overall performance.

.. code-block:: bash

  make config="detailed_timers" # compilation with detailed timers

Additional timers will be shown at the end of the simulation and are also
in ``profile.txt``

----

Create the documentation
^^^^^^^^^^^^^^^^^^^^^^^^^

If you have installed the python module ``sphinx``, you can create the documentation
(which you are currently reading) with:

.. code-block:: bash

   make doc

This creates a local *html* website accessible in your ``build/html/`` folder.

----

.. _installModule:

Install the happi module
^^^^^^^^^^^^^^^^^^^^^^^^

A python module, ``happi``, is provided to view, extract and post-process data from
all the diagnostics.
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
  Instead, we provide the script ``Diagnostics.py`` which is able to find the ``happi``
  module and import it into *python*.

  You may add the following command in your own python script::

    >>> execfile("/path/to/Smilei/scripts/Diagnostics.py")
