GPU offloading
----------------------

To support recent supercomputers, :program:`Smilei` has been ported on GPU (graphical processing units).
Initially built to handle complex graphical output and related to rendering and video games,
GPUs appeared only relatively recently in the HPC ecosystem.
Unlike CPUs, GPUs can do much smaller sets of tasks with massive data throughput.

The most powerful supercomputers are almost all based on GPU acceleration.

----

Supported features
^^^^^^^^^^^^^^^^^^^^

* Both AMD's GPUs and Nvidia's GPUs are supported
* Cartesian geometry in 1D, 2D and in 3D, at order 2 of interpolation
* Diagnostics: ``Field``, ``Probe``, ``Scalar``, ``ParticleBinning``, ``TrackParticles``
* Moving Window
* Boundary conditions for Fields: periodic, reflective and silver-muller are supported (no PML or BM)
* Boundary conditions for Particles: periodic, Reflective, thermal, remove and stop are supported
* Collisions (without ionization or nuclear reactions)

A few key features remain to be implemented:

* AM geometry
* Ionization
* PML boundaries
* Envelope solver
* QED processes

----

Guidelines
^^^^^^^^^^^^^^^^^^^^

Make sure to read the :ref:`compilation documentation<compilationgpu>` and
the :ref:`running documentation<rungpu>`.

On a GPU-equiped cluster, the number of MPI processes (or *tasks*) must be equal to
the number of GPUs, and **each MPI process must be bound to one GPU.**

Splitting your simulation domain works somewhat differently from
CPU simulations. Instead of splitting the simulation box in
many small patches, it is recommended to use **only one patch per GPU**
to obtain the best performance. However, for better memory management,
testing other configurations with a few patches per GPU is encouraged.

