Syntax changes
--------------

In the namelist
^^^^^^^^^^^^^^^


In the post-processing module
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

- In :py:meth:`Field() <Smilei.Field>` and :py:meth:`Probe() <Smilei.Probe>`,
  the argument ``slice`` has been renamed ``average``

- In :py:meth:`ParticleDiagnostic() <Smilei.ParticleDiagnostic>` and :py:meth:`Screen() <Smilei.Screen>`,
  the argument ``slice`` has been renamed ``sum``

- In :py:meth:`Probe() <Smilei.Probe>`, the argument ``average`` requires coordinates
  in code units instead of the indices of the bins.




