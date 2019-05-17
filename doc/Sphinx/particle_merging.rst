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

1. Understand the method
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

1.1 Momentum cell decomposition
""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""

Let us defined some notations first. Momentum norm is called :math:`p` and momentum components
:math:`p_{\alpha}` with :math:`\alpha` equal to x, y or z for each particle.
The number of cells in the direction :math:`\alpha` for the discretization is :math:`N_{\alpha}`.
The discretization step in the direction :math:`\alpha` is called :math:`\Delta_{\alpha}`.

In a position merge cell, step 2 starts by the computation of the minimum :math:`p_{\alpha,min}` and maximum :math:`p_{\alpha,max}` momentum boundaries (also given in :numref:`fig_vranic_particle_merging`).
The boundaries define the momentum space that is then discretized.
The momentum space is divided into momentum cells (of size :math:`\Delta_{\alpha}`) following the discretization (:math:`N_{\alpha}`) given by the user.

In :program:`Smilei`, we use both a spherical discretization geometry for the momentum
discretization and  a Cartesian one as it is the case in :numref:`fig_vranic_particle_merging`.
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

This corresponds to :numref:`fig_spherical_coordinates`.

.. _fig_spherical_coordinates:

.. figure:: _static/spherical_coordinates.png
  :width: 50%

  Spherical coordinates used for the momentum cell discretization.

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

1.2 Merging algorithm for mass macro-particles
""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""

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

  w_t = w_a + w_b \\
  \mathbf{p}_t = w_a \mathbf{p}_a + w_b \mathbf{p}_b \\
  \varepsilon_t = w_a \varepsilon_a + w_b \varepsilon_b

The following energy-momentum relation has to be satisfied for macro-particles a and b:

.. math::
  :label: energy_momentum_relation

  \varepsilon^2 = p^2 + 1

To simplify the problem, Vranic *et al* assume that merged macro-particles
have the same weight :math:`w_a = w_b = w_t / 2`
and same energy :math:`\varepsilon_a = \varepsilon_b = \varepsilon_t / w_t`.

.. _fig_vranic_planar_merging:

.. figure:: _static/vranic_planar_merging.png
  :width: 100%

  View of the plane made by vector :math:`\mathbf{d}` and :math:`\mathbf{p_t}`.
  The corresponding Cartesian frame is given by :math:`(\mathbf{e_1}, \mathbf{e_2}, \mathbf{e_3})`.

As illustrated in :numref:`fig_vranic_planar_merging`, it follows that:

.. math::
  :label: new_momentum_relation

  \mathbf{p}_a +  \mathbf{p}_b = \frac{2 \mathbf{p}_t}{w_t} \\
  \mathbf{p}_{a,\perp} = - \mathbf{p}_{b,\perp} \\
  \mathbf{p}_{a,\parallel} = \mathbf{p}_{b,\parallel} = \mathbf{p_t} / w_t

We all :math:`\omega` the angle betweeb :math:`\mathbf{p_a}` and :math:`\mathbf{p_t}`
so that:

.. math::
  :label: angle_omega

  \cos{\omega} = \frac{\mathbf{p_t}}{w_t \mathbf{p_a}}

We define :math:`\mathbf{d}` the cell direction or location vector.
It represents the location (or the direction in spherical coordinates) of the momentum cell where the macro-particles are located
as shown in :numref:`fig_momentum_cell_vector`.

.. _fig_momentum_cell_vector:

.. figure:: _static/vranic_momentum_cell_vector.png
  :width: 100%

  Momentum cell vector in Cartesian and spherical geometries.

The plane :math:`(\mathbf{e_1},\mathbf{e_2})` is the plane made by the vector :math:`\mathbf{p_t}` and :math:`\mathbf{d}`.
We decide that it contains :math:`\mathbf{p_a}` and :math:`\mathbf{p_b}` so that we have only one possible solution.

Now, it is just necessary to determine :math:`\mathbf{e_1}` and :math:`\mathbf{e_2}` in the momentum frame used by the PIC code.
They are given by the following formula:

.. math::
  :label: planar_coordinates_e1

  \mathbf{e_1} = \mathbf{p_t} / p_t

.. math::
  :label: planar_coordinates_e3

  \mathbf{e_3} & = &  \frac{ \mathbf{d} \times \mathbf{e_1} }{d} \\
               & = & \frac{ 1 }{d.p_t}
   \begin{array}{|l}
      p_{t,z} \cdot d_y - p_{t,y} \cdot d_z \\
      p_{t,x} \cdot d_z - p_{t,z} \cdot d_x \\
      p_{t,y} \cdot d_x - p_{t,x} \cdot d_y
   \end{array}

.. math::
  :label: planar_coordinates_e2

  \mathbf{e_2} & = & \mathbf{e_1} \times \mathbf{e_3} \\
               & = & \frac{1}{p_t^2 . d}
   \begin{array}{|l}
      p_{t,y}^2 .d_x - p_{t,x}(d_y.p_{t,y} + d_z.p_{t,z}) + p_{t,z}^2.d_x \\
      p_{t,z}^2 .d_y - p_{t,y}(d_z.p_{t,z} + d_x.p_{t,x}) + p_{t,x}^2.d_y \\
      p_{t,x}^2 .d_z - p_{t,z}(d_x.p_{t,x} + d_y.p_{t,y}) + p_{t,y}^2.d_z
   \end{array}

Finally, the new macro-particle momentum are:

.. math::
  :label: new_macroparticle_momentum

  \mathbf{p_a} = p_a \left( \cos{\left( \omega \right)} \mathbf{e_1} +  \sin{\left(\omega\right)} \mathbf{e_2} \right) \\
  \mathbf{p_b} = p_b \left( \cos{\left( \omega \right)} \mathbf{e_1} -  \sin{\left(\omega\right)} \mathbf{e_2} \right)

The method is summarized grpahically in :numref:`fig_3d_schematic`.
It has been drawn using Python with Matplotlib.
The Python script in available `here <_static/vranic_geometry.py>`_.

.. _fig_3d_schematic:

.. figure:: _static/vranic_3d_schematics.png
  :width: 100%

  3d view of the different vectors involved in the merging method.

The new macro-particle positions are assigned at the position of one of
the merged macro-particles. We have tested to assign them randomly
or to the first macro-particles of the merged list and we did
not observe any difference.

This algorithm does not work when the total momentum :math:`\mathbf{p}_t` of the macro-particles to be merged
is in the direction of :math:`\mathbf{d}`.
In this case :math:`|| \mathbf{e_3} || = 0` and it is not
possible to determine the system :math:`(\mathbf{e}_1, \mathbf{e}_2, \mathbf{e}_3)`.
In this specific case, the merging is not processed.

1.3 Merging algorithm for macro-photons
""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""

Macro-photons can be merged with the same algorithm.
The only difference is that the momentum norm is equal to the energy :math:`\varepsilon = p`.

When the total momentum :math:`\mathbf{p}_t` is in the direction of :math:`\mathbf{d}`, macro-photons can be merged into a single one contrary to the mass macro-particles since :math:`\varepsilon_t = || \mathbf{p}_t ||`.
This specific situation is implemented in the code.

.. _vranic_implementation:

2. Implementation
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

The Vranic merging method is implemented with the Cartesian
and the Spherical momentum discretization in the source directory ``Merging``.
It is considered as a particle operator and the merging algorithm is managed with a factory (``MergingFactory.h``) as any operator with multiple implementations.
The Cartesian implementation is done in the class ``MergingVranicCartesian`` and the Sphericla one in ``MergingVranicSpherical``.

For both methods, the implemented algorithm is very similar.

    For each cells (in the real space):

    1. Initialization of the momentum cell discretization
    2. Computation of the cell direction vectors (:math:`\mathbf{d}`): this step depends on the discretization and can be efficiently vectorized.
    3. Comutation of the momentum cell indexes for each macro-particle. Efficiently Vectorizable.
    4. Computation of the number of particles per momentum cells.  Not vectorizable because of random memory accesses.
    5. Computation of the cell index of each momentum cell in the sorted array of particles (only the particle indexes are sorted). Not vectorizable.
    6. Sorting of the macro-particles per momentum cells, the cell index previously computed determine where starts each momentum cell. Not vectorizable.

    Then, for each momentum cell:

    1. Division of the macro-particles of the momentum cell in small packs according to the user parameters
    2. Merge of the packs using the previously described Vranic algorithm. Partly vectorized.
    3. Creation of the merged macro-particles at the position of the previous ones
    4. Tag of the macro-particles to be removed

    Then, once the merging finished for a given patch:

    1. Compression of the macro-particle list (remove hole let by removed and tagged particles). By cleaning the particle vector at the end, we limit the computational impact of this step.

2.1 Cartesian momentum Cell discretization
""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""

How to discretize the momentum space is in fact one of the most important point.
The user gives :math:`N_x`, :math:`N_y` and :math:`N_z` via the namelist.
The momentum space boundary corresponds to :math:`p_{\alpha,min}` and :math:`p_{\alpha,max}` with :math:`\alpha` equal to x, y or z.
For this discretization, we force the origin (:math:`p_x = p_y = p_z = 0`) to not be contained in a cell so that there is not in the same cell particles with positive and negative momenta.
The user-defined discretiztion can be slightly adjusted for algorithmic reasons.

    For each momentum component :math:`p_\alpha` with :math:`\alpha` equal to x, y or z:
        If :math:`p_{\alpha,min}` is very close to :math:`p_{\alpha,max}`:
            If :math:`p_{\alpha,min}` and :math:`p_{\alpha,max}` have the same sign:
                Only one cell is used for this component.
                The unique momentum cell is centered around the average particle momentum.
            If :math:`p_{\alpha,min}` and :math:`p_{\alpha,max}` have opposite sign:
                Two cells are used, one for the negative and one for the positive values.
                The discretization is therefore centered in 0.
        Else:
            If :math:`N_\alpha = 1`:
                The unique cell has the size of :math:`p_{\alpha,max} - p_{\alpha,min}`.
            Else if :math:`p_{\alpha,min}` and :math:`p_{\alpha,max}` have the same sign:
                The discretization is classically computed using :math:`N_\alpha`.
            Else if :math:`p_{\alpha,min}` and :math:`p_{\alpha,max}` have opposite sign:
                The discretization is adjusted so that :math:`p_{\alpha} = 0` is at the boundary between 2 consecutive cells. We do it by shifting the discretization and adding an extra cell. At the end, there is an additonal cell than requested (:math:`N_\alpha` = :math:`N_\alpha` + 1).
                

2.2 Spherical momentum Cell discretization
""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""

The user gives :math:`N_r`, :math:`N_\theta` and :math:`N_\phi` via the namelist.
The momentum space boundary corresponds to :math:`p_{r,min}`, :math:`p_{r,max}`, :math:`\theta_{min}`, :math:`\theta_{max}`, :math:`\phi_{min}` and :math:`\phi_{max}`.

    For each momentum component :math:`p_r`, :math:`\theta` and :math:`\phi`:
        If the the minimum boundary is too close to the maximum boundary:
            Only one cell is used for this component.
        Else:
            If :math:`N_\alpha = 1` (here :math:`\alpha` is :math:`r`, :math:`\theta` or :math:`\phi`):

2.3 Solid angle correction
""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""

With the classical spherical discretization, the solid angle that represents the surface crossed by the macro-particles having the same momentum cell direction depends on this direction as shown in :numref:`fig_spherical_discretization` a). In our discretization, the solid angle is larger near :math:`\phi = 0` (equator) and smaller near :math:`\phi = \pi / 2` (poles). Therefore, momentum cells near the equator will potentially have more particles than cells near poles and will undergo more particle merging processes.

.. _fig_spherical_discretization:

.. figure:: _static/spherical_discretization.png
  :width: 100%

  Classical spherical discretization (a) and the spherical discretization with solid angle correction (b). This figure has been generated with the following `Python script <_static/vranic_spherical_discretization.py>`_.

To composate this phenomenon, the discretization (number of cells) in :math:`\theta`, :math:`N_\theta`, is made to depend on :math:`\phi` so that the solid angle is approximatly constant. For this aim, a reference solid angle :math:`\Omega_{ref}` has to be set . It corresponds to the solid angle at the smallest  :math:`|\phi|` value with the :math:`\theta` discretization given by the user in the namelist. For larger :math:`|\phi|` values, the :math:`\theta` discretization :math:`N_\theta` varies to satisfy :math:`\Omega = \sin{(\phi)}\Delta \theta \Delta \phi = \Omega_{ref}`. An example of such a discretization is shown in :numref:`fig_spherical_discretization` b).

2.4 Accumulation effect
""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""

.. _vranic_namelist:

3. Namelist
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

.. _vranic_simulation results:

4. Simulation results
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
