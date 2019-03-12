Vectorization
----------------------

For enhanced performances on most recent CPUs, :program:`Smilei` exploits
efficiently vectorization using refactored and optimized operators.

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

It has been observed that vectorization is more efficient than using scalar
operators when the number of
particles per cell is sufficiently high.
The threshold is evaluated around 10 particles on recent Intel
architectures (Skylake, Knights Landing (KNL), Broadwell).
Vectorization efficiency increases with the number of particles per cell.
Around 256 particles per cell, a speed-up of x2 has been obtained on Intel Skylake
and a speed-up of x3 on Intel KNL using the AVX512 instruction set.
For few particles per cell, scalar implementations are still more efficient
and the ratio is significant.

----

Adaptive vectorization
^^^^^^^^^^^^^^^^^^^^^^^^

Adaptive vectorization consists on locally switch between the scalar and
vectorized operators during the simulation, choosing the most efficient one
in the region of interest.
The concept has been successfully implemented at the lower granularity of the code.
Every given number of time steps, for each
patch, and for each species, the most efficient operator is determined
from the number of particles per cell.

Adaptive vectorization has been validated on large-scale simulations.
One of the case was the simulation of Mildly-relativistic collisionless.
The simulation is illustrated by :numref:`weibel_3d_ne_vecto_it510_fig2`.

.. _weibel_3d_ne_vecto_it510_fig2:

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

The electron density and the patch computational state for the electron species
are shown.
Adaptive vectorization puts the high-density regions rich in
particles in vectorized mode.
Incoming plasma flows, with 8 particles per cell in average, are in scalar mode.
On examined cases, this method allows for speed-ups from x1.3 to x2 regarding only
the macro-particle operators.

This work has been recently submitted for publication
and is avaliable on `ArXiV <https://arxiv.org/abs/1810.03949>`_.
