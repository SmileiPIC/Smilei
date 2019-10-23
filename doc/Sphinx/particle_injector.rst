Particle Injector
================================================================================

Particle injectors enable to inject continuously a flow of fresh macro-particles from the boundaries
into the simulation domain.

1. Understand the method implemented in Smilei
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

In the PIC loop structure, the particle injection is located at the end of the process
after the current and charge projection. It is in the third part of the PIc structure after
the particle sorting and exchange finalization.
            
Injected macro-particles therefore do not contribute to the current and fields of the current iteration
but they are taken into account in the diagnostics.

The method implemented in Smilei in schematically shown in :ref:`fig_particle_injector`.

.. _fig_particle_injector:

.. figure:: _static/figures/particle_injector.png
    :width: 100%

    Description of the particle injection method.
    
For each patch near the boundaries, macro-particles for injectors are first initialized behind domain boundaries.
For the moment, momentum initialization is limited to Maxwellian distributions with a derive.
Macro-particles are pushed using their momentum of a time step :math:`\Delta t` without fields.
Macro-particles that stay behind the domain boundary (due to their momentum direction or because they did not
cross the boundary during this artificial timestep) are not taken into account and removed.
Macro-particles in the domain are kept and considered injected.
They are put in the patch list of macro-particles for the next timestep.

2. Recommendation
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

- Although the code can inject a single species, we recommend to use injectors to inject neutral plasmas.
  This means that positive and negative species should be simultaneously injected at the same boundary.
  To strengthen neutrality, species can be created at the same position.
