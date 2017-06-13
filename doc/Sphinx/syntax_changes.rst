Syntax changes
--------------

In the namelist
^^^^^^^^^^^^^^^


In the post-processing module
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

- The argument ``slice`` replaced by ``average`` or ``sum``, depending on the diagnostic.

- The argument ``stride`` has been replaced by a more complete ``subset``.

- In :py:meth:`Probe() <Smilei.Probe>`, the argument ``average`` requires coordinates
  in code units instead of the indices of the bins.

