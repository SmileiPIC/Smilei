Parallelization basics
----------------------

For high performance, :program:`Smilei` makes complex use of parallel computing,
and it is important to understand the basics of this technology. Parallel simply
means that many processors can run the simulation at the same time, but there is
much more than that.

----

Nodes, cores, processes and threads
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Supercomputers have complex architectures, mainly due to their processors
capability to **work together on the same memory space**. More precisely, *cores*
are grouped in *nodes* (1). All the cores in one node share the same memory space.
This hardware architecture is summarized in :numref:`NodesCoresThreads`.

.. _NodesCoresThreads:

.. figure:: _static/NodesCoresThreads.png
  :width: 11cm
  
  Simplified super-computer architecture.

This same figure shows how the software is structured. *Processes* are the macroscopic
units which are designed to compute over a reserved space in memory (one process
will not handle the memory of another process), and manage several cores.
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
  
  (1) The terminology of *nodes, cores, processes and threads* is not universal. Depending
  on the computer, software (etc.), they can have various meanings.

----

Managing processes and threads
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
Although processes do not share their memory, they must sometimes synchronize their
advance in the execution of the program, or communicate data between each other.
For instance, to calculate the total energy in the simulation, they must communicate
their contribution to the others and compute the sum.
These tasks are accomplished by the Message Passing Interface (MPI) protocol.

At the thread level, the communications do not work in the same manner because threads
already share their data. However, they need synchronization and management to decide
which core handles which thread. This is accomplished by the *OpenMP* protocol.

An illustration of the roles of MPI and OpenMP is provided in :numref:`MPIandOpenMP`

.. _MPIandOpenMP:

.. figure:: _static/MPIandOpenMP.png
  :width: 9cm
  
  MPI handles process-to-process communications, while OpenMP manages threads in a given process.

----

Decomposition of the box
^^^^^^^^^^^^^^^^^^^^^^^^

Traditionally, PIC codes would
split the spatial grid into :math:`N` domains, where :math:`N` is the number
of cores. Each core would manage its own domain on a separate memory space,
and information was communicated between cores using the MPI protocol.
:program:`Smilei` proposes an more efficient approach:
it also decomposes the spatial grid in several domains,
but one core is not directly associated to one domain.

Let us explain this difference in details.
:numref:`PatchDecomposition` gives an example of a grid containing 960 cells.
It is decomposed in :math:`4\times8 = 32` domains, called **patches**.
Each patch has :math:`5\times6` cells.
These patch size is actually reasonable for :program:`Smilei`, whereas
traditional PIC codes would have much larger domains.

The issue is now to decide where these patches will be stored in the memory,
and to choose which cores should do which patches.
Recall that all the cores handled by one process share the same memory:
we will refer to this memory as an *MPI region*.
This means that one process manages one exclusive MPI region.
:numref:`PatchDecomposition` shows an example with the 32 patches split in 5 regions
recognized by their different colors.
Note that these regions are all contiguous, but not necessarily rectangular.

.. _PatchDecomposition:

.. figure:: _static/PatchDecomposition.png
  :width: 10cm
  
  Decomposition of a grid in *patches* and *MPI regions*.

Each MPI region is handled by all the threads of the process. For example, if there are
4 threads in the process that handles the region colored in green, this means the
4 threads will handle 10 patches. The 4 threads will work in parallel, patch by patch,
until all patches are done.

The great advantage of this scheme is that, inside one MPI region, the threads do not
need to wait for their friends to go to the next patch; they can continue working on
the available patches, thus avoiding long waiting times.
This is a form of **local load balancing**.

.. rubric:: Rules

* In each direction :math:`x`, :math:`y`, :math:`z`, the number of patches must be
  a power of 2.
* There must be more patches than threads.


----

Load balancing between MPI regions
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

As we just explained, threads treat the patches in one MPI region asynchronously to
balance their loads carried. Unfortunately, it may not be sufficient.
Indeed, when one MPI region holds much more load than the others, it will take a long
time to compute, while the other processes have already finished and wait for this one.
This can cause large delays.

:program:`Smilei` has an algorithm able to reduce this imbalance by exchanging patches 
from one MPI region to another. A process that has too much load will give patches to
other processes in order to reduce the size of its MPI region. This algorithm is based
on an ordering of the patches by a *Hilbert curve*, as drawn in
:numref:`PatchDecompositionHilbert`. One MPI region contains only patches that contiguously
follow this curve. If this "portion" of the curve has too much load, it will send
some patches to the portions ahead or after, along the same curve. By repeating this
operation every now and then, we ensure that all regions manage an equitable number
of patches. 

.. _PatchDecompositionHilbert:

.. figure:: _static/PatchDecompositionHilbert.png
  :width: 8cm
  
  The shape of the Hilbert curve which determines the patch order.


As the patches can be small, moving a patch from one MPI region to another is
fast: it can fit more easily in the cache, and does not require heavy memory
access.



----

Recommendations
^^^^^^^^^^^^^^^

* **Have as many MPI processes as nodes** in order to optimize the memory distribution.

* On each node, **have as many threads as cores per node**.
  If you have less threads than cores, you will not be using all your cores.
  Use more threads than cores only if hyper-threading is recommended on your architecture.
  
* Use dynamic scheduling for the OpenMP parallelism, by setting the environment variable ``OMP_SCHEDULE``::
    
    export OMP_SCHEDULE=dynamic
    
  This affects only the particles treatment, which will dynamically assign threads.
  Note that fields are always statically assigned to threads.

* **Have small patches**. They can efficiently be as small as 5 cells in each direction.
  This allows good cache use, but also ensures that you have at least as many threads
  as patches, so that they can be treated in parallel.
