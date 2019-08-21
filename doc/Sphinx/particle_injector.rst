Particle Injector
================================================================================

Particle injectors enable to inject continuously a flow of fresh macro-particles from the boundaries
in the simulation domain.

1. Understand the method implemented in Smilei
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

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
            
