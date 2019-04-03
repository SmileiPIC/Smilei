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

The method of M. Vranic basically consists on 3 main steps and is schematically described (in 2D) in Fig. :numref:`fig_vranic_particle_merging`:

1. To Decompose macro-particles in position space into groups so that they share close location. In :program:`Smilei`, macro-particles are sorted by field cells. In the article of M. Vranic *et al.*, the decomposition can be larger than just a cell.

2. Then, to subdivide the macro-particles into sub-groups in the momentum space so that they share close kinteic properties.

3. To merge macro-particles located in the same groups in 2 new macro-particles to respect the charge, energy and momentum conserving laws.

.. _fig_vranic_particle_merging:

.. figure:: _static/vranic_particle_merging.png
  :width: 100%

  Basic description of the Vranic merging method in 2D geometry. In 3D, the idea is strictly the same.

This method has several advantages. It is relatively easy to understand and to implement.
It has a relatively low computational costs and is efficient without
impacting significantly the physical resuls.

.. warning::

  This suppose that the parameters are adequatly tuned.
  Otherwise, the macro-particle merging can affect the final simulation results.

In a position merge cell, step 2 starts by the computation of the minimum and maximum momentum
boundaries (refered to as :math:`p_{min}` and :math:`p_{max}` in :numref:`fig_vranic_particle_merging`).
The boundaries define the momentum space that is then discretized.
The momentum space is divided into momentum cells following a given discretization
given by the user for instance.

In :program:`Smilei`, we use a spherical discretization geometry for the momentum
discretization instead of a Cartesian one as it is the case in :numref:`fig_vranic_particle_merging`.
The momentum space decomposition is basically the same except that the boundaries now concern
the directions :math:`p`, :math:`\theta` and :math:`\phi` in 3D as shown in :numref:`fig_vranic_momentum_discretization`.

.. _fig_vranic_momentum_discretization:

.. figure:: _static/vranic_momentum_discretization.png
  :width: 100%

  2D Cartesian and spherical momentum discretization.

The spherical components are related to the Cartesian momentum components by:

.. math::
  :label: spherical_discretization

  p = \sqrt{ p_x^2 + p_y^2 + p_z^2 }\ ;
  \theta = \arctan{ \left( p_y / p_x \right)}\ ;
  \phi = \arcsin{\left( pz / p \right)}

Since macro-particle momentum components are defined in the Cartesian geometry
by default, considering a spherical discretization induces small additional computation.
However, it makes the merging process more accurate.
Indeed, in the Cartesian discretization, the maximum angle between the momentum
directions of two macro-particle located in the same momentum cell
(i.e. :math:`\theta` and :math:`\phi`) depends on the momentum cell.
For instance, two macro-particles can make an angle up to :math:`\pi / 2` in the cell
adjacent to the origin :math:`p_x = p_y = p_z = 0` whatever the discretization.
In general, this angle diminishes with the distance to the origin.
This issue is therefore negligible for high-energy particles but not
anymore for cold ones.
The spherical geometry ensures that the merging accuracy depends
on the discretization and is similar for all momentum cells.
The overhead induced by the change of geometry is a small fraction of the entire process.

Step 3 starts after the momentum space discretization.
For each momentum cell with more than 4 macro-particles,
the algorithm enables to merge them into 2.
Let us call :math:`\mathrm{M}` the macro-particles in a given momentum cell,
:math:`k` is an index to list each macro-particles of :math:`\mathrm{M}`.
The macro-particle weight is called :math:`w`, the energy :math:`\varepsilon`,
the momentum :math:`\mathbf{p}`.
We start by computing total quantities for the weight :math:`w_t`,
the energy :math:`\varepsilon_t`,
the momentum :math:`\mathbf{p}_t`:

.. math::
  :label: total_quantities

  w_t = \sum_{k \in \mathrm{M}}{w_k}\ ;
  \varepsilon_t = \sum_{k \in \mathrm{M}}{w_k \varepsilon_k}\ ;
  \mathbf{p}_t = \sum_{k \in \mathrm{M}}{w_k \mathbf{p}_k}\ ;

In spherical geometry, the total angles can also be defined:

.. math::
  :label: total_angles

  \theta_t = \sum_{k \in \mathrm{M}}{w_k \theta_k}\ ;
  \phi_t = \sum_{k \in \mathrm{M}}{w_k \phi_k}

To merge all the macro-particles into just one does not allow to locally
conserve weight, energy and momentum. Vranic *et al.* proposes to merge to 2 macro-particles:

.. math::
  :label: merged_particle_relation

  w_t = w_a + w_b \
  \mathbf{p}_t = w_a \mathbf{p}_a + w_b \mathbf{p}_b \
  \varepsilon_t = w_a \varepsilon_a + w_b \varepsilon_b

The following energy-momentum relation has to be satisfied for macro-particles a and b:

.. math::
  :label: energy_momentum_relation

  \varepsilon^2 = p^2 + 1

To simplify the problem, Vranic *et al* assume that merged macro-particles
have the same weight :math:`w_a = w_b = w_t / 2`
and same energy :math:`\varepsilon_a = \varepsilon_b = \varepsilon_t / w_t`.

As illustrated in :numref:`vranic_planar_merging`, it follows table_h_attrs

.. math::
  :label: energy_momentum_relation

  \mathbf{p}_a

Implementation
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Namelist
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Simulation results
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
