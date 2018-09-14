Install
-------

Before installing :program:`Smilei`, you need to install a few dependencies:

* A C++11 compiler, optionally implementing openMP
* an MPI library supporting ``MPI_THREAD_MULTIPLE``
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

#. Download the latest :program:`Smilei` tarball :ref:`here <latestVersion>`.

#. Extract the tarball at the location of your choice.
   Let us assume it is located in your home directory ``~/smilei/``.

#. In a terminal, go to that location and compile:

   .. code-block:: bash

     cd ~/smilei
     make

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
  make print-XXX               # Prints the value of makefile variable XXX
  make env                     # Prints the values of all makefile variables
  make help                    # Gets some help on compilation
  sed -i 's/PICSAR=FALSE/PICSAR=TRUE/g' makefile; make -j4 #To enable calls for PSATD solver from picsar


Each machine may require a specific configuration (environment variables, modules, etc.).
Such instructions may be included, from a file of your choice, via the ``machine`` argument:

.. code-block:: bash

  make machine=my_machine_file

where ``my_machine_file`` is a file, located in ``scripts/CompileTools/machine``, containing
the lines of command to be executed before compilation.

If you successfully write such a file for a common supercomputer, please share it
with developpers so that it can be included in the next release of :program:`Smilei`.

----

Advanced compilation options for profiling/tracing tools
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Several ``make`` options are available in order to compile Smilei ready for
specific performance analysis and debugging tools.

Smilei performance can be analyzed with :
- Scalasca
- Intel Advisor
- Intel Vtune
- Intel Inspector

.. code-block:: bash

  make config="scalasca"             : compilation for scalasca (required the Scalasca profiler)
  make config="advisor"              : compilation for Intel Advisor (required the Intel suite)
  make config="vtune"                : compilation for Intel Vtune (required the Intel suite)
  make config="inspector"            : compilation for Intel Inspector (required the Intel suite)
----

Advanced compilation options for detailed timers
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

The code contains more timers than in the default configuration that can
be activated using a specific option called ``detailed_timers``.
This flag will add the following line to the compilation flags ``-D__DETAILED_TIMERS``
and therefore activate in the code their computation.
These timers are not available in the default compilation
because they are more intrusive and may impact the overall performance in production.

Some of the timers are called inside patches for specific operators
such as the particle pusher. The code first computes the average over the patches
for all MPI domain before computing the final min, max and mean values between
MPI processes.

The final information is therefore an average that does not reflect the
load imbalance between patches of a MPI domain.

the following ``make`` command will activate the additional timers:

.. code-block:: bash

  make config="detailed_timers"      : compilation that activates the detailed timers inside the patch

Additional timers will be shown at the end of the simulation and are also
in ``profile.txt``

----

Compile the documentation
^^^^^^^^^^^^^^^^^^^^^^^^^

If you have installed the python module ``sphinx``, you can compile the documentation
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
