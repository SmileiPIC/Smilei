Units
-----

Just as any PIC code, :program:`Smilei` handles only **dimension-less variables**,
normalized to *reference* quantities.

----

Basic reference quantities
^^^^^^^^^^^^^^^^^^^^^^^^^^

The speed of light, the electron charge and the electron mass provide the basis
of the normalizations in :program:`Smilei`:

* Reference electric charge :math:`Q_r = e` (the electron charge)
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

Normalizing all quantities to these references is convenient for resolving Maxwell's equations,
as it converts them into a dimension-less set of equations.
Indeed, if the let us write :math:`\mathcal{E}=E/E_r`, :math:`\mathcal{B}=B/B_r`,
:math:`\mathcal{J}=j/J_r` and :math:`\mathcal{R}=\rho/N_r/Q_r`, the normalized
electric field, magnetic field, current and density, respectively. We obtain the following
form of Maxwell's equations:

.. math::
  
  \nabla\cdot\mathcal{E} = \mathcal{R}
  \quad\quad
  \nabla\cdot\mathcal{B} = 0
  
  \nabla\times\mathcal{E} = - \partial_t \mathcal{B}
  \quad\quad
  \nabla\times\mathcal{B} = \mathcal{J} + \partial_t \mathcal{E}
  

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
  
  import math
  lambda_ = 2. * math.pi
  cell_length = [0.05 * lambda_]
  sim_length  = [100. * lambda_]


----

Problems requiring explicit units
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Sometimes, :program:`Smilei` may be requested to compute other things than Maxwell's
equations. That is the case, for example, for computing :doc:`collisions <collisions>` or ionization.
In these situations, equations cannot be normalized to dimension-less terms, and
the code must know the value of :math:`\omega_r` in physical units. This requires
defining an :ref:`extra variable in the namelist <wavelength_SI>`.

For instance, ``wavelength_SI = 1e-6`` means that the reference wavelength is one micron,
or equivalently that :math:`L_r = 1\,\mathrm{\mu m} /(2\pi)`. This information will be used only
in some specific parts of the code (collisions, ionization, ...) but not in the main 
PIC algorithms.

.. warning::
  
  The outputs of the code are not converted to SI.
  They are all kept in the reference units listed above.
