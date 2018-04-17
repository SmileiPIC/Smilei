Ionization using a user-defined rate
------------------------------------

:program:`Smilei` now allows for the treatment of ionization considering a user-defined rate.
The Monte-Carlo procedure behind the treatment of ionization in this case closely follows
that developed for :doc:`fied ionization<field_ionization>`.
However, at the moment, only single ionization event per timestep are possible.

The ionization rates are define, for a given ``Species``, as described :ref:`here <Species>`.

----


Benchmarks
^^^^^^^^^^

In what follows, we present two benchmarks of the ionization using a user-defined rate.

The first benchmark considers an initially neutral species that can be potentially ionized twice.
To run this case, a constant and uniform ionization rate is considered that depends only on the particle current charge 
state. For this particular case, we have considered a rate :math:`r_0 = 0.1` (in code units) for ionization from
charge state 0 to 1, and a rate :math:`r_1 = 0.05` (in code units) for ionization from charge state 1 to 2.
The simulation results presented in Fig. :numref:`FigFromRateIonization` (top panel) shows the time evolution of the 
fraction in each possible charge states (:math:`Z=0`, :math:`Z=1` and :math:`Z=2`). 
Super-imposed (dashed lines) are the corresponding theoretical predictions.

The second benchmark features an initially neutral species homogeneously distributed in the simulation box.
The ionization rate is here chosen as a function of the spatial coordinate :math:`x`, 
and reads :math:`r(x) = r_0 \exp(-(x-x_c)^2/2)` with :math:`r_0 = 0.02` the maximum ionization rate and
:math:`x_c=5` the center of the simulation box.
The simulation results presented in Fig. :numref:`FigFromRateIonization` (bottom panel) shows, 
at the end of the simulation :math:`t=20`, the electron number density as a function of space.
Super-imposed (in red) is the corresponding theoretical prediction.

.. _FigFromRateIonization:

.. figure:: _static/userDefinedRate.png
  
  Results of the two benchmarks for the ionization model using user-defined rates as described above.
  




