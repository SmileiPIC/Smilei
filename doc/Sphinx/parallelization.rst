Parallelization
---------------

For high performance, :program:`Smilei` makes complex use of parallel computing,
and it is important to understand the basics of this technology.

----

Nodes, cores, processes and threads
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Supercomputers now have complex architectures, mainly due to their processors
capability to **work together on the same memory space**. More precisely, *cores*
are grouped in *nodes*. All the cores in one node share the same memory space, and two
cores located on two distinct nodes do not have access to the same memory space.
This differs from older computers, where two distinct cores never shared their memory.
This hardware architecture is summarized in :numref:`NodesCoresThreads`.

.. _NodesCoresThreads:

.. figure:: _static/NodesCoresThreads.png
  :width: 14cm
  
  Simplified super-computer architecture.

The same figure shows how the software is structured. *Processes* are the macroscopic
units which are designed to compute over a reserved space in memory (one process
will not handle the memory of another process), and handle a various number of cores.
To treat these cores, a process must provide several *threads*. A thread is basically the
sequence of instructions from the program, which must be executed by one core.
However, a thread is not uniquely associated to one core: a core can run two threads,
one after the other.

The association between the software *threads* and the hardware *cores* can be more
complicated. :numref:`NodeWith2Processes` shows an example where two processes share the
same node. In this case, we illustrate the memory of this node as split in two parts because
the two processes cannot access to the same memory.

.. _NodeWith2Processes:

.. figure:: _static/NodeWith2Processes.png
  :width: 6cm
  
  An example where two processes share the same node.

.. warning::
  
  The terminology of *nodes, cores, processes and threads* is not universal. Depending
  on the computer, software (etc.), they can have various meanings.

----

Managing processes and threads
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Although processes do not share their memory, they must sometime communicate data or
synchronize their advance in the execution of the program. For instance, they may have to
wait for other processes so that they all are at the same point of the simulation.
Another example: to calculate the total energy in the simulation, they must communicate
their contribution to the others and compute the sum. These tasks are accomplished by
the *Message Passing Interface* (MPI) protocol.

At the thread level, the communications do not work in the same manner because threads
already share their data. However, they need synchronization and management to decide
which core handles which thread. This is accomplished by the *OpenMP* protocol.

An illustration of the roles of MPI and OpenMP is provided in :numref:`MPIandOpenMP`

.. _MPIandOpenMP:

.. figure:: _static/MPIandOpenMP.png
  :width: 10cm
  
  MPI handles process-to-process communications, while OpenMP manages threads in a given process.

----

Decomposition of the box
^^^^^^^^^^^^^^^^^^^^^^^^

:program:`Smilei` proposes an innovative approach
to decompose the simulation between all these elements. Traditionally, PIC codes would
split the spatial grid into :math:`N` domains, where :math:`N` is the number
of cores. Each core would handle its own domain on a separate memory space,
and information was communicated between domains using the MPI protocol.
:program:`Smilei` also decomposes the spatial grid in several
domains, but one core is not directly associated to one domain.

Let us explain this difference in details.
:numref:`PatchDecomposition` gives an example of a grid containing 960 cells.
It is decomposed in :math:`4\times8 = 32` domains, called **patches**.
These patches can be very small: down to :math:`5\times5` cells. 

The issue is now to decide where these patches will be stored in the memory.
Recall that all the cores handled by one process share the same memory:
we will refer to this memory as an *MPI region*.
A clever algorithm assigns several patches to each MPI region.
:numref:`PatchDecomposition` shows an example with the 32 patches split in 5 regions
recognized by their different colors.
Note that these regions are not necessarily rectangular.

.. _PatchDecomposition:

.. figure:: _static/PatchDecomposition.png
  :width: 10cm
  
  Decomposition of a grid in *patches* and *MPI regions*.

Each MPI region is handled by all the threads of the process. For example, if there are
4 threads in the process that handles the region colored in green, this means the
4 threads will handle 10 patches. The 4 threads will successively work on the next
available patch until all patches are done.

.. rubric:: Advantages

* Inside one MPI region, the threads do not need to wait for their friends to go to the
  next patch; they can continue working on the next patch, thus avoiding long waiting
  times. This is a form of **local load balancing**.
* The patches are regularly exchanged between MPI regions in order to uniformize the load
  carried by each. This is a form of **global load balancing**.
* As the patches can be small, moving a patch from one MPI region to another is
  fast: it can fit more easily in the cache, and does not require heavy memory
  access.


.. rubric:: Rules

* In each direction :math:`x`, :math:`y`, :math:`z`, the number of patches must be
  a power of 2.
* There must be more patches than threads.


.. rubric:: Recommendations

* **Have as many MPI processes as nodes** in order to optimize the memory sharing.
* On each node, **have as many threads as cores per node**.
  If you have less threads than cores, you will not be using all your cores.
  Use more threads than cores only if hyper-threading is recommended on your architecture.
* Use dynamic scheduling for the OpenMP parallelism, by setting the environment variable ``OMP_SCHEDULE``::
    
    export OMP_SCHEDULE=dynamic
    
  This variable affects only the particles treatment, which will dynamically assign threads.
  Note that fields are always statically assigned to threads.
* **Have small patches**. They can efficiently be as small as 5 cells in each direction.
  This allows good cache use, but also ensures that you have at least as many threads
  as patches, so that they can be treated in parallel.

.. rst-class:: inprogress
  
  In progress ...
