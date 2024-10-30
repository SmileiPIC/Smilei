Particle Injector
================================================================================

Particle injectors provide a continuous flow of fresh
macro-particles from the boundaries into the simulation domain.

----

The method implemented in Smilei
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

In the PIC loop structure, the particle injection occurs
after current projection on the grid, particle sorting and synchronizations.
Injected macro-particles therefore do not contribute to the current and fields
of the current iteration but they are taken into account in the diagnostics.

.. _fig_particle_injector:

.. figure:: /_static/figures/particle_injector.png
    :width: 100%


----

Recommendation
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

* Although a single species may be injected, we recommend to inject both
  positively and negatively charged species at the same time to ensure
  a neutral plasma. To strengthen neutrality, species may be created at
  the same positions.

* If the particle momentum is drawn from a Maxwellian, we recommend to use a random
  positionning instead of the regular one.
  Regular positionning may induce numerical effects such as loss of charge and spurious field near the boundary.
  The reason is explained in the following figure.
  The regular positionning works when injecting a drifting cold plasma with a drift velocity
  sufficiently high to let the particles entering the simulation domain.

.. _fig_particle_injector_regular_random:

.. figure:: /_static/figures/particle_injector_regular_random.png
    :width: 100%

----

Implementation
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

The particle injector algorithm is coded in the file
``Patch/VectorPatch.cpp`` in the function ``injectParticlesFromBoundaries``.

The class ``ParticleInjector`` manages the injector's parameters and properties,
while new macro-particles are initialized using the class ``ParticleCreator``.
