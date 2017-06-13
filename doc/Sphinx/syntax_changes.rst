Syntax changes
--------------

In the namelist
^^^^^^^^^^^^^^^

- The ``DiagParticles`` diagnostic has been renamed :ref:`DiagParticleBinning <DiagParticleBinning>`.


In the post-processing module
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

- The ``ParticleDiagnostic()`` method has been renamed :py:meth:`ParticleBinning() <Smilei.ParticleBinning>`.

- The argument ``slice`` replaced by ``average`` or ``sum``, depending on the diagnostic.

- The argument ``stride`` has been replaced by a more complete ``subset``.

- In :py:meth:`Probe() <Smilei.Probe>`, the argument ``average`` requires coordinates
  in code units instead of the indices of the bins.

