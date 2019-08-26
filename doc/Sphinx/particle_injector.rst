Particle Injector
================================================================================

Particle injectors enable to inject continuously a flow of fresh macro-particles from the boundaries
into the simulation domain.

1. Understand the method implemented in Smilei
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

In the PIC loop structure, the particle injection is located at the end of the process
after the projection and exchanges.

.. code-block:: ReST

    for each patch:
        for each species:
            
            Interpolate
            Pusher
            
        for each species:
        
            Projection
            
    Maxwell
            
    for each patch:
        for each species:
        
            Exchange finalization
            Cleaning
            Particle sorting
            **Particle injection**
            
    Diagnostics
            
Injected macro-particles therefore do not contribute to the current and fields of the current iteration
but they are taken into account in the diagnostics.

The method implemente in Smilei in shcematicllay shown in :ref:`fig_particle_injector`

.. _fig_particle_injector:

.. figure:: _static/figures/particle_injector.png
    :width: 100%

    Description of the particle injection method.
