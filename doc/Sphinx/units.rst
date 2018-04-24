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
  
  Counter-intuitively, the reference density :math:`N_r` is **not equal** to :math:`L_r^{-3}`.

Normalizing all quantities to these references is convenient for resolving Maxwell's equations,
as it converts them into a dimension-less set of equations:

.. math::
  
  \mathbf{\nabla}\cdot\mathbf{E} = \rho
  \quad\quad
  \nabla\cdot\mathbf{B} = 0

  \nabla\times\mathbf{E} = - \partial_t \mathbf{B}
  \quad\quad
  \nabla\times\mathbf{B} = \mathbf{j} + \partial_t \mathbf{E}

where :math:`\mathbf{E}`, :math:`\mathbf{B}`, :math:`\mathbf{j}` and :math:`\mathbf{\rho}`
are the electric field, magnetic field, current density and charge density, normalized to
:math:`E_r`, :math:`B_r`, :math:`J_r` and :math:`Q_r N_r`, respectively. Note that the
temporal and spatial derivatives are also normalized to :math:`T_r` and :math:`L_r`, respectively.

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

.. _Weights:

Macro-particle weights
^^^^^^^^^^^^^^^^^^^^^^

Macro-particles are assigned a *statistical weight* which represents
their contribution to the plasma distribution function. 
In :program:`Smilei`, this weight is defined at the beginning of the simulation
for each particle and is never modified afterwards. Its definition reads:

.. math::
  
  \textrm{macro-particle weight} = \frac
      {\textrm{species density in the cell}}
      {\textrm{number of particles of this species in the cell}}

As a consequence, the sum of all weights of the particles in one cell is equal to
the density of the species in this cell, in units of :math:`N_r`.

The charge carried by a macro-particle (in Coulomb) can therefore be retrieved by the inverse operation:

.. math::

 \textrm{macro-particle charge} = \tilde{w}\tilde{q}V_{\rm cell}N_r


:math:`\tilde{w}` and :math:`\tilde{q}` are the normalized weight and charge of the particle.
:math:`V_{\rm cell}` is the volume of a single cell of the simulation.
This last equation can be used to derive charge carried by track particles or to initialize weights in the case of a species which positions are initialized via a numpy array (see :ref:`Species`).


