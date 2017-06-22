Syntax changes
--------------

In the namelist
^^^^^^^^^^^^^^^

- Parameter ``maxwell_sol`` renamed :py:data:`maxwell_solver`.

- Parameters ``bc_em_type_*`` unified in a new parameter :py:data:`EM_boundary_conditions`.

- Parameter ``referenceAngularFrequency_SI`` renamed :py:data:`reference_angular_frequency_SI`.

- Parameter ``initPosition_type`` renamed :py:data:`position_initialization`

- Parameter ``initMomentum_type`` renamed :py:data:`momentum_initialization`

- Parameters ``bc_part_type_****`` unified in a new parameter :py:data:`boundary_conditions`. 

- Particle boundary conditions do not accept ``"none"`` anymore. Use ``"periodic"`` instead.

- Parameter ``thermT`` renamed :py:data:`thermal_boundary_temperature`

- Parameter ``thermVelocity`` renamed :py:data:`thermal_boundary_velocity`

- Parameter ``isTest`` renamed :py:data:`is_test`

- Parameter ``boxSide`` renamed :py:data:`box_side`

- Parameter ``polarizationPhi`` renamed :py:data:`polarization_phi`

- Parameter ``dump_file_sequence`` renamed :py:data:`keep_n_dumps`

- Diagnostic ``DiagParticles`` renamed :ref:`DiagParticleBinning <DiagParticleBinning>`.

- Removed parameters ``track_******``, and created the :ref:`DiagTrackParticles <DiagTrackParticles>` block.


In the post-processing module
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

- ``ParticleDiagnostic()`` method  renamed :py:meth:`ParticleBinning() <Smilei.ParticleBinning>`.

- Argument ``slice`` replaced by ``average`` or ``sum``, depending on the diagnostic.

- Argument ``stride`` replaced by a more complete ``subset``.

- In :py:meth:`Probe() <Smilei.Probe>`, the argument ``average`` requires coordinates
  in code units instead of the indices of the bins.

