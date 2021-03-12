Single-domain multiple decompositions
-------------------------------------

In Single Domain Multiple Decompositions (SDMD) the paralelization slightly differs from a standard simulation.
When this technique is used, numerical objects called regions and dedicated to field operations are created.
There are as many regions as MPI processes and each process owns exactly one region.
A single region extent therefore covers the equivalent of many patches allowing communicationless operations over larger area.
Operations requiring many ghost cells, like large stencil filters or spectral solvers, can therefore be much more efficiently executed.
It also makes up for the synchronization overhead induced by the use of very small patches during times of the simulation where the load is well balanced.
The benefits of this option are illustrated `in this paper <https://hal.archives-ouvertes.fr/hal-02973139>`_.

Particles are still handled as in a standard decomposition. 
Fields are simply communicated from regions to patch in order to interpolate their values at the particles locations.
Inversely, currents and densities are communicated from patches to regions after deposition.

The regions extents are always Cartesian and as homogeneous in size as possible
The association between MPI processes and regions are dynamically reconsidered during the simulation if dynamic load balancing is used.
The goal is to maximize the overlap between the region and the collection of patches owned by the MPI process at all times in order to reduce communication overhead.
