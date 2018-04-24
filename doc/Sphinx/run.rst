.. _run:

Run
---

Before you launch :program:`Smilei`, :doc:`write a namelist file<namelist>`
containing all the information of your simulation (grid shape, particles, lasers, diagnostics, etc.).

You can also start from an example provided in the ``benchmarks`` directory.

----

The ``smilei`` executable
^^^^^^^^^^^^^^^^^^^^^^^^^

`Compiling Smilei <compile>`_ creates an executable file ``smilei`` in the source directory.

.. code-block:: bash
  
  smilei arg1 arg2 arg3 ...

The command-line arguments ``arg1``, ``arg2``, ``arg3`` (etc.) can be:

* the path to a namelist
* any *python* instruction that you want to execute during the namelist reading.

The simplest example, to run your namelist ``my_namelist.py``, is

.. code-block:: bash
  
  ./smilei  my_namelist.py

You may also add an additional instruction to be appended at the end of the namelist:

.. code-block:: bash
  
  ./smilei  my_namelist.py  "Main.print_every=10"

Note that, in addition, you will generally use the ``mpirun`` or ``mpiexec`` command
to run :program:`Smilei` on several MPI processes:

.. code-block:: bash
  
  mpirun -n 4 ./smilei  my_namelist.py  "Main.print_every=10"

If you want to run several openMP threads per MPI processes, you usually have to set
the following environment variable to the desired number of threads before running
``mpirun``:

.. code-block:: bash
  
  export OMP_NUM_THREADS=4

When running :program:`Smilei`, the output log will remind you how many MPI processes and openMP threads
your simulation is using.

----

Running in *test mode*
^^^^^^^^^^^^^^^^^^^^^^

A second executable ``smilei_test`` is available (after the usual compilation)
to run in the *test mode*:

.. code-block:: bash
  
  ./smilei_test my_namelist.py

This *test mode* does the same initialization as the normal mode,
except it only loads the first patch of the full simulation. After initialization,
the test mode exits so that the PIC loop is *not* computed.

This mode may be used to check the consistency of the namelist, and to make sure
simple errors will not occur. It does not check all possible errors, but it runs fast.

Running in **test mode requires to run on 1 MPI process only**. However, it is possible
to indicate what is the partition of MPI processes and OpenMP threads intended for the
actual simulation. For instance, to test your namelist that is intended to run on 1024 MPI
processes, each hosting 12 OpenMP threads, use the following syntax:

.. code-block:: bash
  
  ./smilei_test 1024 12 my_namelist.py

----

Directory management
^^^^^^^^^^^^^^^^^^^^

Let us assume you have written your namelist ``my_namelist.py``, and that you placed it
inside your home directory. Also, we assume that the ``Smilei`` directory is also there,
so that the ``smilei`` executable is located in ``~/Smilei/``.

Knowing that :program:`Smilei` generally writes out all the results in the current directory,
it is recommended to create a new directory to store these results. For instance:

.. code-block:: bash
  
  $ mkdir ~/my_simulation                     # New directory to store results
  $ cp ~/my_namelist.py ~/my_simulation       # Copies the namelist there
  $ cd ~/my_simulation                        # Goes there
  $ mpirun -n 4 ~/Smilei/smilei my_namelist   # Run with 4 processors

----

Using the provided script
^^^^^^^^^^^^^^^^^^^^^^^^^

For simple cases such as the previous one, use the script ``smilei.sh``, provided in
the `Smilei` directory. You only have to run

.. code-block:: bash
  
  $ ./smilei.sh 4 my_namelist.py

where the number 4 says that the code will run 4 MPI processes. A directory with all
the results will automatically be created next to your namelist.

----

Running on large clusters
^^^^^^^^^^^^^^^^^^^^^^^^^

We do not provide instructions to run on super-computers yet. Please refer to your
administrators.


----

Debugging
^^^^^^^^^

In case of problems, the code can be compiled with additional debugging flags (usual ``-g`` and ``-O0``) and internal 
checks by compiling it with 

.. code-block:: bash
  
    make config=debug

Compiling the whole code with this command will make it very slow to run. 
But to check only a particular file for errors, first compile the code with ``make``, then
modify the file, and recompile in debug mode.

In debug mode, these C++ macros are activated:

* ``DEBUG("some text" [<< other streamable])``
* ``HEREIAM("some text" [<< other streamable])``


----

Known issues
^^^^^^^^^^^^

* OpenMPI ``2.*`` often causes unstable behavior in Smilei.
  For instance, with ``openmpi 2.1``, the `vader` protocol seems to interfere with Smilei's
  memory management and comunications. We therefore recommend to disable this
  protocol when running ``mpirun``, as follows:

  .. code-block:: bash
  
    $ mpirun --mca btl ^vader -n 4 ~/Smilei/smilei my_namelist   # Disable vader

