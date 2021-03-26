.. rst-class:: experimental

Single-domain multiple decompositions
-------------------------------------

The standard decomposition of the simulated domain consists in splitting
the grid in rectangular "patches" that contain both fields and particles
(see :doc:`parallelization`). A new technique has been developed for Smilei,
"Single-Domain Multiple Decompositions" (SDMD), where two decompositions are
made:

* Particles are still contained in the same small, rectangular patches.
* Fields are contained in large, rectangular "regions".

Each MPI process owns many patches but only (and exactly) one region.
A single region extent therefore covers the equivalent of many patches allowing
communication-less operations over a larger area.

Fields are simply communicated from regions to patches in order to interpolate their
values at the particles locations. Inversely, currents and densities are communicated
from patches to regions in order to solve Maxwell's equations.

The region extents are always Cartesian and their size is as homogeneous as possible.
The association between MPI processes and regions is dynamically reconsidered
during the simulation if dynamic load balancing is used.
The goal is to maximize the overlap between the region and the collection of
patches owned by the MPI process at all times in order to reduce communication overhead.

.. rubric:: Advantages

* The synchronization overhead induced by the use of very small patches is reduced.
* Operations requiring many ghost cells (like large-stencil filters or spectral solvers)
  can be executed much more efficiently.

The benefits of SDMD are illustrated
`in this paper <https://hal.archives-ouvertes.fr/hal-02973139>`_.
