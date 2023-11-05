GPU offloading
----------------------

In order to be able to be used on most recent supercomputers :program:`Smilei` has been ported on GPUs.
Although a few key features remain to be implemented such as AM geometry, ionization and PML, the fundamentals of the code are currently accelerated on GPUs
if GPUs are available.

Initially built to handle complex graphical output and related to rendering and video games, GPUs appeared only relatively recently in the HPC ecosystem.
Unlike CPUs, GPUs can do much smaller sets of tasks with massive data throughput.

Currently 7 of the 10 most powerful supercomputers are based on GPU acceleration and everything points to that trend of adding accelerators to continue: The announced exaflopic supercomputers will include GPUs.


* Guideline using GPU: in general it is better to use one patch per GPU to obtain the best performance out of the GPU  

* Current supported features:

  * Both AMD's GPUs and Nvidia's GPUs are supported
  * Cartesian geometry in 2D and in 3D
  * Diagnostics: field, probe, scalar, bin, particle tracking 
  * Moving Window.

* Notable features currently missing: AM geometry, ionization, PML, envelop
