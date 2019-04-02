Particle Merging
================================================================================

The ability to merge macro-marticles can speed-up the code efficiency
and reduce the memory footprint in some specific simulation senarii:

* When macro-particles accumulate in a fraction of the simulation domain
  hence strongly worsening the load imbalance (ex: Weibel collison shocks,
  laser wakefield electron acceleration).
* When macro-particles are generated in a large quantity due to some
  additonal physical mechanisms (ionization, macro-photon emission, QED pair production...)
* When macro-particles travel in large quantities outside interesting physical regions.

Available implemented methods:

* M. Vranic merging method (`M. Vranic et al., CPC, 191 65-73 (2015) <https://doi.org/10.1016/j.cpc.2015.01.020>`_)

--------------------------------------------------------------------------------

.. _ref_vranic_method:

The merging method of M. Vranic
--------------------------------------------------------------------------------

.. _ref_understand_vranic_method:

Understand the method
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

The method of M. Vranic basically consists on 3 main steps and is schematically described (in 2D) in Fig. :numref:`vranic_merging_method`:

1. To Decompose macro-particles in position space into groups so that they share close location. In :program:`Smilei`, macro-particles are sorted by field cells. In the article of M. Vranic *et al.*, the decomposition can be larger than just a cell.

2. Then, to subdivide the macro-particles into sub-groups in the momentum space so that they share close kinteic properties.

3. To merge macro-particles located in the same groups in 2 new macro-particles to respect the charge, energy and momentum conserving laws.

.. figure:: _static/vranic_particle_merging.png
  :width: 90%

  Basic description of the Vranic merging method in 2D geometry. In 3D, the idea is strictly the same.

This method has several advantages. It is relatively easy to understand and to implement.
It has a relatively low computational costs and is efficient without
impacting significantly the physical resuls.

.. warning::

  This suppose that the parameters are adequatly tuned.
  Otherwise, the macro-particle merging can affect the final simulation results.

In a position merge cell, step 2 starts by the computation of the minimum and maximum momentum
boundaries (refered to as :math:`p_{min}` and :math:`p_{max}` in Fig. :numref:`vranic_merging_method`).
The boundaries define the momentum space that is then discretized.
The momentum space is divided into momentum cells following a given discretization
given by the user for instance.

In :program:`Smilei`, we use a spherical discretization geometry for the momentum
discretization instead of a Cartesian one as shown in Fig. :numref:`vranic_merging_method`.

.. figure:: _static/momentum_discretization.png
  :width: 90%

  Basic description of the Vranic merging method in 2D geometry. In 3D, the idea is strictly the same.


Implementation
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Namelist
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Simulation results
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
