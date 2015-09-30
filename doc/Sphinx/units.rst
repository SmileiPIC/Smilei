Units
-----

Just as any PIC code, :program:`Smilei` handles only **dimension-less variables**,
normalized to *reference* quantities. The most basic ones are: 

* the reference electric charge :math:`Q_r = e` (the electron charge)
* the reference mass :math:`M_r = m_e` (the electron mass)
* the reference velocity :math:`V_r = c` (the speed of light)
* the reference energy :math:`K_r = m_e c^2`
* the reference momentum :math:`P_r = m_e c`

Even with these normalizations, :program:`Smilei` **does not know the scale of the problem**:
it lacks a reference distance, or equivalently, a reference time.
Instead of choosing a physical constant (for example, the electron radius), the scale of
the problem is not decided *a priori*, and the user is free to scale the result of the
simulation to any value.
In fact, quantities are proportional an *unknown* reference frequency
:math:`\omega_r`, which can be scaled by the user *a posteriori*.
Usually, :math:`\omega_r` will be an important frequency of the problem,
for example, a laser frequency, if there is one.

From this reference frequency :math:`\omega_r`, we define:

* a reference time :math:`T_r = 1/\omega_r`
* a reference length :math:`L_r = c/\omega_r` 
* a reference electric field :math:`E_r = m_e c \omega_r / e`
* a reference magnetic field :math:`B_r = m_e \omega_r / e`
* a reference charge density :math:`N_r = \varepsilon_0 m_e \omega_r^2 /e^2`
* a reference current :math:`J_r = c\, e\, N_r`

Normalizing all quantities to these references is convenient for resolving Maxwell's equations,
as it converts them into a dimension-less set of equations.

In the :doc:`namelist <namelist>`, the user must provide all parameters in units of :math:`Q_r`, :math:`M_r`,
:math:`V_r`, :math:`K_r`, :math:`P_r`, :math:`T_r`, :math:`L_r`, :math:`E_r`, :math:`B_r`,
:math:`N_r` or :math:`J_r`.

Sometimes, :program:`Smilei` may be requested to compute other things than Maxwell's
equations. That is the case, for example, for computing :doc:`collisions <collisions>` or ionization.
In this situation, :program:`Smilei` must know the value of :math:`\omega_r` in physical units,
and defining an :ref:`extra variable in the namelist <wavelength_SI>` is required.

