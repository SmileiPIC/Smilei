.. _multiphotonBreitWheelerPage:

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
is given by the Ritus formulae

.. math::
  :label: BWEnergyDistribution

  \frac{d^2N_{BW}}{d \chi_{\pm} dt} = \frac{\alpha_f m_e^2 c^4}{\pi \sqrt{3} \hbar \varepsilon_\gamma \chi_\gamma}
  \int_{x}^{+\infty}{\sqrt{s} K_{1/3} \left( \frac{2}{3} s^{3/2} \right) ds - \left( 2 - \chi_\gamma x^{3/2} \right) K_{2/3} \left( \frac{2}{3} x^{3/2} \right) }

where :math:`x = \left( \chi_\gamma / (\chi_{-} \chi_{+}) \right)^{2/3}`.
The parameters :math:`\chi_{-}` and :math:`\chi_{+}` are the respective Lorentz
invariant of the electron and the positron after pair creation.
Furthermore, one has :math:`\chi_- = \chi_\gamma - \chi_+` meaning that :math:`\chi_-`
and :math:`\chi_+` can be interchanged.

The total production rate of pairs can be written

.. math::
  :label: BWproductionRate

  \frac{dN_{BW}}{dt} = \frac{\alpha_f m_e^2 c^4}{ \hbar \varepsilon_\gamma} \chi_\gamma T \left( \chi_\gamma \right)

where

.. math::
  :label: BWTfunction

  T \left( \chi_\gamma \right) = \frac{1}{\pi \sqrt{3} \chi_\gamma^2 }
  \int_{0}^{+\infty}{\int_{x}^{+\infty}{\sqrt{s} K_{1/3} \left( \frac{2}{3} s^{3/2}
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

The right subplot in :numref:`synchrotron_pairs_dNdt` gives the probability
for a photon to decay into a pair as a function of the energy given to the electron
(using approximation :math:`\chi_\gamma / \chi_- = \gamma_\gamma / \gamma_-`)
for a field of amplitude :math:`E = 500 m_e \omega c / e`.
It can also be seen as the pair creation energy distribution.
The distribution is symmetric with respect to :math:`\chi_- / \chi_\gamma = 1/2`.
Below :math:`\chi_\gamma = 10`, The maximum probability corresponds to
equal electron-positron energies :math:`\chi_- = \chi_+ = \chi_\gamma / 2`.
Above this threshold, the energy dispersion increases with :math:`\chi_\gamma`.

.. _synchrotron_pairs_dNdt:

.. figure:: _static/synchrotron_pairs_dNdt.png
  :width: 18cm

  (left) - Normalized total pair production distribution given by Eq. :eq:`BWproductionRate`.
  (right) - Normalized pair creation :math:`\chi` distribution given by Eq. :eq:`BWEnergyDistribution`.


--------------------------------------------------------------------------------

.. _BWStochasticSchemeSection:

Stochastic scheme
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

The Multiphoton Breit-Wheeler is treated with a Monte-Carlo process similar
to the nonlinear inverse Compton Scattering
(see :ref:`the radiation reaction page <radiationReactionPage>`).
It is close to what has been done in
[Duclous2011]_, [Lobet2013]_, [Lobet2015]_.

The first preliminary step consists
on introducing the notion of macro-photon. Macro-photons are simply the equivalent of
macro-particles (see :ref:`the macro-particle section <QuasiParticlesSection>`)
extended to photons.
There are defined by a charge and a mass equal to 0. The momentum is substituted
by the photon momentum :math:`\mathbf{p}_\gamma = \hbar \mathbf{k}` where
:math:`\mathbf{k}` is the wave vector.
The momentum contains the photon energy so that
:math:`\mathbf{p}_\gamma = \gamma_\gamma m_e \mathbf{c}`.
The definition of the photon Lorentz factor is therefore also slightly different
than particles.

1. An incremental optical depth :math:`\tau`, initially set to 0, is assigned to the macro-photon.
Decay into pairs occurs when it reaches the final optical depth :math:`\tau_f`
sampled from :math:`\tau_f = -\log{(\xi)}` where :math:`\xi` is a random number in :math:`\left]0,1\right]`.

2. The optical depth :math:`\tau` evolves according to the photon quantum parameter
following:

.. math::
  :label: mBW_MCDtauDt

  \frac{d\tau}{dt} = \frac{dN_{BW}}{dt}\left( \chi_\gamma \right)

that is also the production rate of pairs
(integration of Eq. :eq:`BWEnergyDistribution`).

3. The emitted electron's quantum parameter :math:`\chi_-` is computed by
inverting the cumulative distribution function:

.. math::
  :label: mBW_CumulativeDistr

  P(\chi_-,\chi_\gamma) = \frac{\displaystyle{\int_0^{\chi_-}{
  \frac{d^2N_{BW}}{d \chi dt} d\chi}}}{\displaystyle{\int_0^{\chi_\gamma}{\frac{d^2N_{BW}}{d \chi dt} d\chi}}}

The inversion of  :math:`P(\chi_-,\chi_\gamma)=\xi'` is done after drawing
a second random number
:math:`\xi' \in \left[ 0,1\right]` to find :math:`\chi_-`.
The positron quantum parameter is :math:`\chi_+ = \chi_\gamma - \chi_-`.

4. The energy of the emitted electron is then computed:
:math:`\varepsilon_- = mc^2 \gamma_- = mc^2 \left[ 1 + \left(\gamma_\gamma - 2\right) \chi_- / \chi_\gamma \right]`.
If :math:`\gamma_\gamma < 2`, the pair creation is not possible since the photon
energy is below the rest mass of the particles.

5. The photon momentum is then updated.
Propagation direction is the same as for the photon. Pairs are created at the
same position as for the photon. The weight is conserved. It is possible to
create more than a macro-electron or a macro-positron in order to improve
the phenomenon statistics. In this case, the weight of each macro-particle is
the photon weight divided by the number of emissions.

--------------------------------------------------------------------------------

Implementation
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

C++ classes for the multiphoton Breit-Wheeler process are located
in the directory ``src/MultiphotonBreitWheeler``.
In :program:`Smilei`, the multiphoton Breit-Wheeler process is not incorporated
in the photon pusher in order to preserve vector performance of the latter one.

Description of the files:

* Class ``MultiphotonBreitWheelerTables``: this class contains the methods to generate the tables,
  to output them, to read them and to broadcast them among MPI tasks.
  It also contains methods to get values from the tables for the Monte-Carlo process.
* Class ``MultiphotonBreitWheeler``: this class contains the methods to
  perform the Breit-Wheeler Monte-Carlo process described in :ref:`the previous section <BWStochasticSchemeSection>`).
* Class ``MultiphotonBreitWheelerFactory``: this class is supposed to
  manage the different Breit-Wheeler algorithms.
  For the moment, only one model is implemented.

Formula :eq:`BWTfunction` and :eq:`mBW_CumulativeDistr` are tabulated
at the beginning of the simulation because of the cost of their computation
for each photon.
The parameters such as the table ranges and discretization can be
given in the :ref:`MultiphotonBreitWheeler <MultiphotonBreitWheeler>` namelist section.
Once generated, the table can be written on the disk and reloaded for a next run.
Small tables coded in hdf5 are provided in the repository in the folder
databases with the name: `multiphoton_Breit_Wheeler_tables.h5`.

If the multiphoton Breit-Wheeler is activated for a photon species, the factory
will initialize the instance ``Multiphoton_Breit_Wheeler_process`` of
the class ``MultiphotonBreitWheeler``
declared in the corresponding ``species`` (see ``species.cpp``).

The multiphoton Breit-Wheeler Monte-Carlo process is performed in the method ``dynamics`` of ``species``.
It is called after the particle field interpolation (field gathering),
after ionization and radiation reaction and before the particle pusher.
At this stage, the new particles are stored in a temporary buffer called ``new_pair``.
This is an array of two instances of ``Particles``.
It is declared in ``Multiphoton_Breit_Wheeler_process``.
Particles are imported in the main species particle arrays
(``particles`` object in ``species``) only after the current deposition
and before the boundary conditions using the method ``importParticles``
of the class ``Particles``.

--------------------------------------------------------------------------------

Benchmarks
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Synchrotron, 2D
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

Counter-propagating plane wave, 1D
""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""

.. _counter_pair_smilei:

.. figure:: _static/counter_pair_smilei.png
  :width: 18cm

  (left) - Energy balance of the simulation.
  (middle) - Final energy spectrum of the electrons (blue), positrons (orange), and photons (green).
  (right) - Time evolution of the number of macro-electrons (blue),
  macro-positrons (orange) and macro-photons (green) in the simulation.
