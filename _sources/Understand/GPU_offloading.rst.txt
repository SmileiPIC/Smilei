GPU offloading
----------------------

To support recent supercomputers, :program:`Smilei` has been ported on GPU (graphical processing units).
Initially built to handle complex graphical output and related to rendering and video games,
GPUs appeared only relatively recently in the HPC ecosystem.
Unlike CPUs, GPUs can do much smaller sets of tasks with massive data throughput.

Currently 7 of the 10 most powerful supercomputers are based on GPU acceleration and
the trend seems confirmed at least in the near future: 
the announced exaflopic supercomputers will include GPUs.

* Guideline: in general it is recommended to use one patch per GPU
  to obtain the best performance. However, for better memory management,
  testing a few patches per GPU is encouraged.

* Currently supported features:

  * Both AMD's GPUs and Nvidia's GPUs are supported
  * Cartesian geometry in 1D, 2D and in 3D , for order 2
  * Diagnostics: Field, Probes, Scalar, ParticleBinning, TrackParticles
  * Moving Window
  * Boundary conditions for Fields: Periodic, reflective and silver-muller are supported (no PML or BM)
  * Boundary conditions for Particles: Periodic, Reflective, thermal, remove and stop are supported

* A few key features remain to be implemented (AM geometry, ionization, PML, envelope,
  additional physics), but the fundamentals of the code are ported.

* Short term roadmap We are currently working on porting on GPU the following features: AM geometry (order 2) and collisions
