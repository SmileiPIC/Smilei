Parallelization
---------------

Supercomputers now have complex architectures, mainly due to the computing units
(CU) able to work together on the same memory space. More precisely, CUs
are grouped in *nodes*. All the CUs in one node share the same memory space, and two
CUs located on two distinct nodes do not have access to the same memory space.
This differs from older supercomputers, where two distinct CUs never shared their memory.

To take advantage of this new technology, :program:`Smilei` proposes a new approach
to decompose the simulation between all the CUs. Traditionally, PIC codes would
split the spatial grid into :math:`N_{CU}` domains, where :math:`N_{CU}` is the number
of CUs. Each CU would handle its own domain, and information was communicated between domains
using the MPI protocol. :program:`Smilei` also decomposes the spatial grid in several
domains, but one CU is not directly associated to one domain.

Let us explain this difference in details.
:numref:`PatchDecomposition` gives an example of a grid containing 960 cells.
It is decomposed in :math:`4\times8 = 32` sub-domains, called **patches**.
These patches can be very small: down to :math:`5\times5` cells. 

The issue is now to decide where these patches will be stored in the memory.
Recall that all the CUs in one *node* share the same memory: we will refer to this
memory as an *MPI region*. A clever algorithm assigns several patches to each MPI
region. :numref:`PatchDecomposition` shows an example with the 32 patches split
in 5 regions recognized by their different colors.
Note that these regions are not necessarily rectangular.

.. _PatchDecomposition:

.. figure:: _static/PatchDecomposition.png
  :width: 10cm
  
  Decomposition of a grid in *patches* and *MPI regions*.

Each MPI region is handled by all the CUs of the node. For example, if there are
4 CUs in the node that handles the region colored in green, this means the 4 CUs will
handle 10 patches. The 4 CUs will successively work on the next available patch until
all patches are done.

To handle patches simultaneously in a given node, :program:`Smilei` uses the 
*OpenMP* protocol. This ensure that all the CUs in one node can run together on the
same memory.
  
.. rubric:: Advantages

* Inside one MPI region, the CUs do not need to wait for their friends to go to the
  next patch; they can continue working on the next patch, thus avoiding long waiting
  times. This is a form of **local load balancing**.
* The patches are regularly
  exchanged between MPI regions in order to uniformize the load carried by each.
  This is a form of **global load balancing**.
* As the patches can be small, moving a patch from one MPI region to another is
  fast: it can fit more easily in the cache, and does not require heavy memory
  access.


.. rubric:: Rules

* In each direction :math:`x`, :math:`y`, :math:`z`, the number of patches must be
  a power of 2.
* There must be more patches than CUs.


.. rubric:: Recommendations

* **Have as many MPI processes as nodes** in order to optimize the memory sharing.
* In each node, **have as many OpenMP threads as CUs**. If you have less threads than CUs,
  you will not be using all your CUs. If you have more threads than CUs (overthreading),
  some threads will be treated sequentially.
* **Have small patches**. They can efficiently be as small as 5 cells in each direction.
  This allows good cache use, but also ensures that you have at least as many threads
  as patches, so that they can be treated in parallel.

.. rst-class:: inprogress
  
  In progress ...
