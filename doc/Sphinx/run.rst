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

Note that, in addition, you will generally use the ``mpiexec`` or ``mpirun`` command
to run :program:`Smilei` on several processors:

.. code-block:: bash
  
  mpiexec -np 4 ./smilei  my_namelist.py  "Main.print_every=10"

----

Running in *test mode*
^^^^^^^^^^^^^^^^^^^^^^

A ``-T`` option may be provided to the ``smilei`` executable to run in the *test mode*.
It must be *the first argument* provided:

.. code-block:: bash
  
  ./smilei -T my_namelist.py

This *test mode* does the same initialization as the normal mode,
except it only loads the first patch of the full simulation. After initialization,
the test mode exits so that the PIC loop is *not* computed.

This test mode may be used to check the consistency of the namelist, and to make sure
simple errors will not occur. It does not check all possible errors, but it runs fast.

Running in **test mode requires to run on 1 MPI process only**. However, it is possible
to indicate what is the partition of MPI processes and OpenMP threads intended for the
simulation. For instance, to test your namelist that is intended to run on 1024 MPI
processes, each hosting 12 OpenMP threads, use the following syntax:

.. code-block:: bash
  
  ./smilei -T1024x12 my_namelist.py

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
  $ mpiexec -np 4 ~/Smilei/smilei my_namelist # Run with 4 processors

----

Using the provided script
^^^^^^^^^^^^^^^^^^^^^^^^^

For simple cases such as the previous one, use the script ``smilei.sh``, provided in
the `Smilei` directory. You only have to run

.. code-block:: bash
  
  $ ./smilei.sh 4 my_namelist.py

where the number 4 says that the code will run on 4 processors. A directory will all
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
But to check only a particular file for errors, first compile the code with `make`, then
modify the file, and recompile in debug mode.

In debug mode, these C++ macros are activated:

* ``DEBUG("some text" [<< other streamable])``
* ``HEREIAM("some text" [<< other streamable])``


----

Reporting bugs
^^^^^^^^^^^^^^

To report bugs, please create an issue on the `github page <https://github.com/SmileiPIC/Smilei/issues/new>`_ .
