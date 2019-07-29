Vectorization
----------------------

For enhanced performances on most recent CPUs, :program:`Smilei` exploits
efficiently vectorization using refactored and optimized operators.

Vectorization optimizations are published in [Beck2019]_.

----

Notion of Single Instruction Multiple Data (SIMD) Vectorization
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Single Instruction Multiple Data (SIMD) vectorization consists on performing on
a contiguous set of data, usually called vector, the same operation(s)
in a single instruction.
On modern Computational Processing Units (CPU), vector registers have a length 512 kb
that corresponds to 8 double precision floats (on Intel Skylake processors for
instance and future ARM architecture).
Each processing unit can perform a Fused Multiply Add instruction (FMA) that
combines an addition and a multiplication.
If-conditions can be handled using mask registers.
Modern SIMD vectorization is described in :numref:`simd_fig`.

.. _simd_fig:

.. figure:: _static/SIMD.png
    :width: 90%
    :align: center

    Single Instruction Multiple Data (SIMD) vectorization

On SIMD CPUs, an application has to use SIMD vectorization to reach the maximum
of the core computational peak performance. A scalar code without FMA
uses less than 7% of the core computational power.
This affirmation can nonetheless be mitigated on Intel Skylake processors that
adapt their frequency on the used vectorization instruction set.

----

SIMD vectorization of the particle operators
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Optimization efforts have been recently done to vectorize efficiently the
particle operators of :program:`Smilei`.

A new sorting method has been first implemented in order to then make
the particle operator vectorization easier.
This method, referred to as cycle sort, minimizes the number of data movements
by performing successive permutation.

The most expensive operators and most difficult to vectorize are the current projection
(deposition) and the field interpolation (gathering) steps where
there is an interpolation between the grids and the macro-particles.
These two steps have been vectorized taking advantage of the cycle sort.

----

Vectorization Performance
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Vectorization is not always the most efficient choice.
It depends on the number of macro-particles per cell.
To demonstrate this, we have evaluated in [Beck2019]_ the performance with a series of tests on different architectures: Intel Cascade
Lake, Intel Skylake, Intel Knights Landing, Intel Haswell, Intel Broadwell.
The Cascade Lake processor is not in the original study and has been added after.
We have used the 3D homogeneous Maxwellian benchmark.
The number of macro-particles per cell is varied from 1 to 512.
This study has been focused on the particle operators (interpolator, pusher, projector, sorting) and discards the
computational costs of the Maxwell solver and of the communications between processes.
Each run has been performed on a single node with both the scalar and the vectorized operators..
Since the number of cores varies from an architecture
to another, the runs were conducted so that the load per core
(i.e. OpenMP thread) is constant.
The number of patches per core also remains the same for all cores throughout the whole simulation since the imbalance
in this configuration is never high enough to trigger patch exchanges.
The patch size is kept constant at 8 × 8 × 8 cells.
The total number of patches for each architecture is determined so that each core has 8 patches to handle.
The numerical parameters are given in :numref:`vecto_numerical_parameters`.

.. _vecto_numerical_parameters:

+-------------------------------------+-------------------------------------------------------+-------------------+---------------------------+
| Cluster                             | Architecture                                          | Number of patches | Configuration             |
+=====================================+=======================================================+===================+===========================+
| Jean Zay, IDRIS, France             | 2 x Cascade Lake (Intel® Xeon® Gold 6248, 20 cores)   | 5 x 8 x 8         | Intel 19, IntelMPI 19     |
+-------------------------------------+-------------------------------------------------------+-------------------+---------------------------+
| Irene Joliot-Curie, TGCC, France    | 2 x skylake (Intel® Skylake 8168, 24 cores)           | 6 x 8 x 8         | Intel 18, IntelMPI 18     |
+-------------------------------------+-------------------------------------------------------+-------------------+---------------------------+
| Frioul, Cines, France               | 2 x Knights Landing (Intel® Xeon® Phi 7250, 68 cores) | 8 x 8 x 8         | Intel 18, IntelMPI 18     |
+-------------------------------------+-------------------------------------------------------+-------------------+---------------------------+
| Tornado, LPP, France                | 2 x Broadwell (Intel® Xeon® E5-2697 v4, 16 cores)     | 4 x 8 x 8         | Intel 17, openMPI 1.6.5   |
+-------------------------------------+-------------------------------------------------------+-------------------+---------------------------+
| Jureca, Juelich, Germany            | 2 x Haswell (Intel® Xeon® E5-2680 v3, 12 cores)       | 3 x 8 x 8         | Intel 18, IntelMPI 18     |
+-------------------------------------+-------------------------------------------------------+-------------------+---------------------------+

The results of the simulation tests (shape factor of order 2) for both scalar and vectorized versions are
shown in :numref:`vecto_particle_times_o2_all`.
Contrary to the scalar mode, the vectorized operators efficiency depends strongly on the number of particles per cell.
It shows improved efficiency, compared to the scalar mode, above a certain number of particles per cell denoted *inversion point*.

.. _vecto_particle_times_o2_all:

.. figure:: _static/vecto_particle_times_o2_all.png
  :width: 100%

  Particle computational cost as a function of the number of particles per cell. Vectorized
  operators are compared to their scalar versions on various cluster
  architectures. Note that the Skylake compilations accepts both AVX512 and AVX2
  instruction sets.

The lower performances of the vectorized operators at low particles per cell can be easily understood:

1. The complexity of vectorized algorithms is higher than their scalar counter-parts.
#. New schemes with additional loops and local buffers induced an overhead that is onmy compensated when the number of particles is large enough.
#. SIMD instructions are not efficient if not fulfilled
#. SIMD instructions operate at a lower clock frequency than scalar ones on recent architectures

The location of the inversion point of the speed-ups brought by vectorization depends on the architecture.
The performance results are summarized in :numref:`vecto_performance_results`.

.. _vecto_performance_results:

+-------------------------------------+-------------------------------------------------------+------------------------+
| Architecture (Cluster)              | Inversion point (particles per cell)                  | Vectorization speed-up |
+=====================================+=======================================================+========================+
| Cascade lake (Jean Zay)             | 8 particles per cell                                  | x2                     |
+-------------------------------------+-------------------------------------------------------+------------------------+
| Skylake (Irene Joliot-Curie)        | 10 particles per cell (most advanced instruction set) | x2.1                   |
+-------------------------------------+-------------------------------------------------------+------------------------+
| KNL (Frioul)                        | 12 particles per cell                                 | x2.8                   |
+-------------------------------------+-------------------------------------------------------+------------------------+
| Broadwell (LLR)                     | 10 particles per cell                                 | x1.9                   |
+-------------------------------------+-------------------------------------------------------+------------------------+
| Haswell (Jureca)                    | 10 particles per cell                                 | x1.9                   |
+-------------------------------------+-------------------------------------------------------+------------------------+

Vectorization efficiency increases with the number of particles per cell above the inversion point.
It tends to stabilize far from the inversion point above 256 particles per cell.


----

Adaptive vectorization
^^^^^^^^^^^^^^^^^^^^^^^^

Adaptive vectorization consists on switching localy between scalar and
vectorized operators during the simulation, choosing the most efficient one
in the region of interest.
The concept has been successfully implemented at the lower granularity of the code.
Every given number of time steps, for each
patch, and for each species, the most efficient set of operator is determined
from the number of particles per cell.
The concept is schematically described in :numref:`fig_vecto_domain_decomposition`.

.. _fig_vecto_domain_decomposition:

.. figure:: _static/vecto_domain_decomposition.png
  :width: 100%

  Description of the adaptive vectorization withn the multi-stage domain decomposition.
  Patches with many macro-particles per cell are faster in with vectorized operators whereas with few macro-particles per cell, scalar operators are more efficient.


Large-scale simulations
^^^^^^^^^^^^^^^^^^^^^^^^

Adaptive vectorization has been validated on large-scale simulations with
different benchmarks.

One of the case was the simulation of Mildly-relativistic collisionless shock.
The effect of the adaptive vectorization mode is illustrated by :numref:`fig_weibel_3d_ne_vecto_it510`.
The electron density is shown in the volume rendering of the top.
The volume rendering at the bottom shows and patch computational state for the electron species.

.. _fig_weibel_3d_ne_vecto_it510:

.. figure:: _static/Weibel_3d_ne_vecto_it510.jpg
    :width: 90%
    :align: center
    :target: https://youtu.be/-ENUekyE_A4

    Mildly-relativistic collisionless shock: On the top, volume rendering of the normalized
    electron density :math:`n_e /n_c` (:math:`n_c` the critical density) at
    time :math:`t = 34 \omega^{-1}` (:math:`\omega` the laser frequency) after the beginning of the collision.
    On the bottom, patches in vectorized
    mode for the electron species at the same time.
    An animated version of these can be viewed by clicking on this image.


Thanks to the adaptive vectorization, high-density regions that contains many macro-particles per cell corresponds to the patches in vectorized mode.
Incoming plasma flows, with 8 particles per cell in average, are in scalar mode.
How the simulation configures dynamically in time the scalar and vectorized regions is demonstrated in the following video:

.. _video_weibel_3d_ne_vecto_it510:

.. raw:: html

  <video style="display:block; margin: 0 auto; width: 100%;" controls src="http://www.maisondelasimulation.fr/projects/Smilei/uploads/videos/weibel_interp.mp4" width="100%">
  </video>

On examined cases, this method allows for speed-ups from x1.3 to x2 regarding only
the macro-particle operators.

----

References
^^^^^^^^^^

.. [Beck2019] `A. Beck et al., ArXiV 1810.03949 (2019) <https://arxiv.org/abs/1810.03949>`_
