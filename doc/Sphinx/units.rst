Units
-----

Like many PIC codes, :program:`Smilei` handles only **dimension-less variables**,
normalized to *reference* quantities.

----

Basic reference quantities
^^^^^^^^^^^^^^^^^^^^^^^^^^

The speed of light, the elementary charge and the electron mass provide the basis
of the normalizations in :program:`Smilei`:

* Reference electric charge :math:`Q_r = e` (the elementary charge)
* Reference mass :math:`M_r = m_e` (the electron mass)
* Reference velocity :math:`V_r = c` (the speed of light)

We can derive from these:

* a reference energy :math:`K_r = m_e c^2`
* a reference momentum :math:`P_r = m_e c`

Even with these normalizations, :program:`Smilei` **does not know the scale of the problem**:
it lacks a reference distance, or equivalently, a reference time.

----

Arbitrary reference quantities
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Instead of choosing a physical constant (for example, the electron radius) as a reference,
the scale of the problem is not decided *a priori*, and the user is free to scale the result
of the simulation to any value.
In fact, quantities are proportional an *unknown* reference frequency
:math:`\omega_r`, which can be scaled by the user *a posteriori*.

Usually, :math:`\omega_r` will be an important frequency of the problem.
For example, if there is a laser, it could be the laser frequency. 
Or it could be the electron plasma frequency.

From this reference frequency :math:`\omega_r`, we define:

* a reference time :math:`T_r = 1/\omega_r`
* a reference length :math:`L_r = c/\omega_r` 
* a reference electric field :math:`E_r = m_e c \omega_r / e`
* a reference magnetic field :math:`B_r = m_e \omega_r / e`
* a reference particle density :math:`N_r = \varepsilon_0 m_e \omega_r^2 /e^2`
* a reference current :math:`J_r = c\, e\, N_r`

.. warning::
  
  :math:`1/N_r` is a volume, but counter-intuitively, it is **not equal** to :math:`L_r^{3}`.

Normalizing all quantities to these references is convenient for resolving Maxwell's equations,
and the charges equation of motion, as it converts them into a dimension-less set of equations:

.. math::

  \mathbf{\nabla}\cdot\mathbf{E} = \rho
  \quad\quad
  \nabla\cdot\mathbf{B} & = 0 \\

  \nabla\times\mathbf{E} = - \partial_t \mathbf{B}
  \quad\quad
  \nabla\times\mathbf{B} = & \; \mathbf{j} + \partial_t \mathbf{E} 

.. math::

  \partial_t \mathbf{p} = Z \mathbf{E} + Z \mathbf{v}\times\mathbf{B}
  
where :math:`\mathbf{E}`, :math:`\mathbf{B}`, :math:`\mathbf{j}` and :math:`\mathbf{\rho}`
are the electric field, magnetic field, current density and charge density, normalized to
:math:`E_r`, :math:`B_r`, :math:`J_r` and :math:`Q_r N_r`, respectively. :math:`Z` and
:math:`\mathbf p` are a particle's charge and momentum, normalized to :math:`Q_r` and 
:math:`P_r`, respectively. Note that the temporal and spatial derivatives are also
normalized to :math:`T_r` and :math:`L_r`, respectively.


----

Tips for the namelist
^^^^^^^^^^^^^^^^^^^^^

In the :doc:`namelist <namelist>`, the user must provide all parameters in units of :math:`Q_r`,
:math:`M_r`, :math:`V_r`, :math:`K_r`, :math:`P_r`, :math:`T_r`, :math:`L_r`, :math:`E_r`,
:math:`B_r`, :math:`N_r` or :math:`J_r`.

This may be cumbersome if you know your input data in other units.
However, the namelist is actually a *python* code that can compute conversions easily.

For example, let us assume that you know your problem size in units of the wavelength.
Knowing that the reference wavelength is :math:`2\pi L_r`, you can multiply all your
lengths by :math:`2\pi`::
  
  from math import pi
  wavelength = 2. * pi
  cell_length = [0.05 * wavelength]
  grid_length  = [100. * wavelength]


----

Problems requiring explicit units
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Sometimes, :program:`Smilei` may be requested to compute other things than Maxwell's
equations. That is the case, for example, for computing :doc:`collisions <collisions>` or ionization.
In these situations, equations cannot be normalized to dimension-less terms, and
the code must know the value of :math:`\omega_r` in physical units. This requires
defining an :ref:`extra parameter in the namelist <reference_angular_frequency_SI>`.

For instance, ``reference_angular_frequency_SI = 2.*pi*3e8/1e-6`` means that
:math:`L_r = 1\,\mathrm{\mu m} /(2\pi)`.
This information will be used only in some specific parts of the code (collisions, ionization, ...)
but not in the main PIC algorithms.

.. warning::
  
  The outputs of the code are not converted to SI.
  They are all kept in the reference units listed above.

----

.. _integrated_quantities:

Quantities integrated over the grid
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Special care must be taken when considering local quantities that are spatially
integrated.

.. rubric:: 1. The spatially-integrated kinetic energy density

The particle kinetic energy density is naturally in units of :math:`K_r N_r`.
Integrating over space give different results depending on the simulation dimension.
In 1D, this space a length, with units :math:`L_r`; in 2D, it is a surface, with units
:math:`L_r^2`; and in 3D, it is a volume, with units :math:`L_r^3`.
Overall, the integrated energy has the units :math:`K_r N_r L_r^D`
where :math:`D` is the simulation dimension. Note that we could expect
to obtain, in 3D, an energy with units :math:`K_r`, but counter-intuitively
it has the units :math:`K_r N_r L_r^3`.

These kinetic energies appear, for instance, in the :ref:`DiagScalar` as
``Ukin`` (and associated quantities).

.. rubric:: 2. The spatially-integrated electromagnetic energy density

The electromagnetic energy density has the units :math:`E_r^2/\varepsilon_0 = K_r N_r`.
Consequently, the spatially-integrated electromagnetic energy density has
the units :math:`K_r N_r L_r^D`; the same as the integrated kinetic energy density above.

These electromagnetic energies appear, for instance, in the :ref:`DiagScalar` as
``Uelm`` (and associated quantities).

.. rubric:: 3. The space- & time-integrated Poynting flux

The Poynting flux has the units :math:`E_r B_r / \mu_0 = V_r K_r N_r`.
Consequently, the flux integrated over a boundary, and over time, has the units
:math:`V_r K_r N_r L_r^{D-1} T_r = K_r N_r L_r^D`, which is the same as the
integrated energy densities above.

This integrated Poynting flux appears, for instance, in the :ref:`DiagScalar` as
``Uelm_bnd``, ``PoyXmin``, ``PoyXminInst`` (and associated quantities).


----

.. _Weights:

Macro-particle weights
^^^^^^^^^^^^^^^^^^^^^^

Macro-particles are assigned a *statistical weight* which represents
their contribution to the plasma distribution function. 
In :program:`Smilei`, this weight is defined at the beginning of the simulation
for each particle and is never modified afterwards. Its definition reads:

.. math::
  
  \textrm{macro-particle weight} = \frac
      {\textrm{species density} \times \textrm{cell hypervolume}}
      {\textrm{number of macro-particles in cell}}

As the density is in units of :math:`N_r` and the cell hypervolume in
units of :math:`L_r^D` (where :math:`D` is the simulation dimension),
then the units of weights is :math:`N_r L_r^D`.

This definition of weights ensures that they do not depend on the
cell hypervolume, i.e. they can be reused in another simulation, as long as
:math:`D`, :math:`L_r` and :math:`N_r` are unchanged.


