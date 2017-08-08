Multiphoton Breit-Wheeler pair creation
--------------------------------------------------------------------------------

The multiphoton Breit-Wheeler (:math:`\gamma + n\omega \rightarrow e^- + e^+`),
also referred to as
the nonlinear Breit-Wheeler, corresponds to the decay of a
high-energy photon into a pair of electron-positron
when interacting with a strong electromagnetic field.

In the vacuum, the electromagnetic field becomes nonlinear from the Schwinger
electric field :math:`E_s = 1.3 \times 10^{18}\ \mathrm{V/m}` corresponding
to an intensity of :math:`10^{29}\ \mathrm{Wcm^{-2}}` for
:math:`\lambda = 1\ \mu m`. In such a field, spontaneous apparitions of electron-positron pairs from
the nonlinear vacuum are possible. If this field is not reachable in the laboratory
frame, it will be very close in the boosted frame of highly-accelerated electrons. At
:math:`10^{24}\ \mathrm{Wcm^{-2}}`, a Lorentz factor of :math:`\gamma \sim 10^5`
is required to reach the Schwinger limit.
This is the reason why quantum electrodynamics effects and in particular strong-
field pair generation is accessible with the extreme-intensity lasers (multipetawatt lasers).

As for electrons or positrons, the strength of QED effects depends
on the photon quantum parameter:

.. math::
  :label: photonQuantumParameter

  \chi_\gamma = \frac{\gamma_\gamma}{E_s} \sqrt{ \left( \mathbf{E}_\perp
  + \mathbf{c} \times \mathbf{B} \right)^2 }

where
:math:`\gamma_\gamma = \varepsilon_\gamma / m_e c^2` is the photon normalized energy,
:math:`m_e` the electron mass,
:math:`c` the speed of light in vacuum,
:math:`\mathbf{E}_\perp` is the electric field orthogonal to
the propagation direction of the photon,
:math:`\mathbf{B}` the magnetic field.

--------------------------------------------------------------------------------

Physical model
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

The energy distribution of the production rate of pairs by a hard photon
is given by

.. math::
  :label: BWEnergyDistribution

  \frac{dN_{BW}}{d \chi_{\pm} dt} = \frac{\alpha_f m_e^2 c^4}{\pi \sqrt{3} \hbar \varepsilon_\gamma \chi_\gamma}
  \int_{0}^{+\infty}{\sqrt{s} K_{1/3} \left( \frac{2}{3} s^{3/2} \right) ds - \left( 2 - \chi_\gamma x^{3/2} \right) K_{2/3} \left( \frac{2}{3} x^{3/2} \right) }

where :math:`x = \left( \chi_\gamma / (\chi_{-} \chi_{+}) \right)^{2/3}`.
The parameters :math:`\chi_{-}` and :math:`\chi_{+}` are the respective Lorentz
invariant of the electron and the positron after pair creation.
Furthermore, one has :math:`chi_- = \chi_\gamma - \chi_+` meaning that :math:`\chi_-`
and :math:`\chi_+` can be interchanged.

The total production rate of pairs can be written

.. math::
  :label: BWproductionRate

  \frac{dN_{BW}}{dt} = \frac{\alpha_f m_e^2 c^4}{ \hbar \varepsilon_\gamma} \chi_\gamma T \left( \chi_\gamma \right)

where

.. math::
  :label: BWTfunction

  T \left( \chi_\gamma \right) = \frac{1}{\pi \sqrt{3} \chi_\gamma^2 }
  \int_{0}^{+\infty}{\int_{0}^{+\infty}{\sqrt{s} K_{1/3} \left( \frac{2}{3} s^{3/2}
  \right) ds - \left( 2 - \chi_\gamma x^{3/2} \right) K_{2/3} \left( \frac{2}{3} x^{3/2} \right) }} d\chi_-

A photon of energy :math:`\varepsilon_\gamma` traveling in a constant electric field :math:`E` has a Lorentz
parameter equal to :math:`\chi_\gamma = \varepsilon_\gamma E / (E_s m_e c^2)`.

We consider the case where photon interact in a constant uniform electric field.
For a field of amplitude :math:`E = 500 m_e \omega c / e`, the energy
distribution and the production rate of pair creation as a function of :math:`\chi_\gamma` are
plotted in :numref:`synchrotron_pairs_dNdt`. It shows that the total production
rate of electron-positron pairs rises rapidly to reach a peak around
:math:`\chi_\gamma = 10` with almost a pair generated per femtosecond.
Under :math:`\chi_\gamma = 0.1`, the production rate is very weak with
less than a pair after 100 picoseconds of interaction.
Above :math:`\chi_\gamma = 10`, the production decreases slowly with
:math:`\chi_\gamma`.

The right plot in :numref:`synchrotron_pairs_dNdt` gives the probability
for a photon to decay into a pair as a function of the energy given to the electron
(using approximation :math:`\chi_\gamma / \chi_- = \gamma_\gamma / \gamma_-`).
It can also be seen as the pair creation energy distribution.
Below :math:`\chi_\gamma = 10`, the distribution is symmetric with respect
to :math:`\chi_- / \chi_\gamma = 1/2`.

.. _synchrotron_pairs_dNdt:

.. figure:: _static/synchrotron_pairs_dNdt.png
  :width: 18cm

  (left) - Normalized total pair production distribution given by Eq. :eq:`BWproductionRate`.
  (right) - Normalized pair creation :math:`\chi` distribution given by Eq. :eq:`BWEnergyDistribution`.

Stochastic scheme
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Implementation
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Benchmarks
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Counter-propagating plane wave, 1D
""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""

Synchrotron 2D
""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""

In this configuration, a mono-energetic bunch of photons is initialized
in a constant uniform strong magnetic field.
The photons decay into pairs via the multiphoton Breit-Wheeler progressively.
In this particular case, the generated electrons and positrons do not radiate
in order to capture the emission energy spectrum.
Two cases are simulated with different
initial quantum parameters:

* Case 1: :math:`\chi_{\gamma,0} = 1`, :math:`B = 270`, :math:`\gamma_{\gamma,0} = 1500`
* Case 2: :math:`\chi_{\gamma,0} = 20`, :math:`B = 600`, :math:`\gamma_{\gamma,0} = 8125`

.. _synchrotron_pairs_energy_spectra_chi1:

.. figure:: _static/synchrotron_pairs_energy_spectra_chi1.png
  :width: 18cm

  (left) - Electron energy spectrum at the end of the run.
  (middle) - Positron energy spectrum at the end of the run.
  (right) - Time evolution of the photon (green), electron (blue)
  and positron (orange)
  normalized energy :math:`U / U_{tot}`.

.. _synchrotron_pairs_energy_spectra_chi20:

.. figure:: _static/synchrotron_pairs_energy_spectra_chi20.png
  :width: 18cm

  (left) - Electron energy spectrum at the end of the run.
  (middle) - Positron energy spectrum at the end of the run.
  (right) - Time evolution of the photon (green), electron (blue)
  and positron (orange)
  normalized energy :math:`U / U_{tot}`.
