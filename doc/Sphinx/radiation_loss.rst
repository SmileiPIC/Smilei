Synchrotron-like radiation loss
--------------------------------------------------------------------------------

I. Introduction
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

High-energy particles traveling in a strong electromagnetic field loss energy by
radiation. Depending on the field strength and the particle energy, radiation
losses occur in different regimes and can be seen as smooth or brutal with
diffusive and stochastic consequences.
This phenomenon is usually referred to as synchrotron-like radiation emission
(in reference to the emission process occurring in synchrotron facilities
with a constant magnetic field)
or nonlinear inverse Compton scattering (arbitrary electromagnetic field).

In extremely intense laser fields such as attainable intensities with future
multi-petawatt facilities above :math:`10^{21}\ \mathrm{Wcm^{-2}}`, high-energy
radiation emission strongly influence the
dynamics of charged particles and the overall energy balance of laser-plasma
interaction.

Different approaches have been implemented in :program:`Smilei` as summarized
in tab :numref:`radiationRegimes` to deal with different regimes of emission.
These regimes can be characterized via the quantum Lorentz invariant parameter
:math:`\chi_\pm` which is an indicator of how strong is the radiation emission
process.

.. math::
  :label: particleQuantumParameter

  \chi_{\pm} = \frac{\gamma_{\pm}}{E_s} \left| \left( \beta \cdot \mathbf{E}
  \right)^2 - \left( \mathbf{E} + \mathbf{v} \times \mathbf{B} \right)^2
  \right|^{1/2}

where :math:`E_s = m^2 c^3 / \hbar e \simeq 1.3 10^{18}\ \mathrm{V/m}` is
the Schwinger field, :math:`e` is the electron charge,
:math:`m` is the electron mass, :math:`c` the speed of light in vacuum,
:math:`\hbar` the reduced Planck constant. :math:`\mathbf{E} = (E_x, E_y, E_z)`
and :math:`\mathbf{B} = (B_x, B_y, B_z)` are respectively the electric and
the magnetic field. :math:`\gamma = \varepsilon_\pm / m c^2` is the particle
Lorentz factor and also the normalized particle energy. :math:`\beta = v/c` is
the normalized particle velocity.

.. _radiationRegimes:

+-------------------------------------+--------------------------+------------------------------------------------+---------------------------+
| Regime                              | :math:`\chi` value       | Description                                    | Models                    |
+=====================================+==========================+================================================+===========================+
| Classical radiation emission        | :math:`\chi \sim 10^{-3}`| :math:`\varepsilon_\gamma  \ll \varepsilon_\pm`| Landau-Lifshitz           |
|                                     |                          | , radiated energy overestimated for            |                           |
|                                     |                          | :math:`\chi > 10^{-2}`                         |                           |
+-------------------------------------+--------------------------+------------------------------------------------+---------------------------+
| Semi-classical radiation emission   | :math:`\chi \sim 10^{-2}`| :math:`\varepsilon_\gamma  \ll \varepsilon_\pm`| Corrected Landau-Lifshitz |
|                                     |                          | , no stochastic effects                        |                           |
+-------------------------------------+--------------------------+------------------------------------------------+---------------------------+
| Weak quantum regime                 | :math:`\chi \sim 10^{-1}`| :math:`\varepsilon_\gamma  < \varepsilon_\pm`, | Stochastic model of       |
|                                     |                          | :math:`\varepsilon_\gamma / mc^2  \gg 1`       | Niel `et al` / Monte-Carlo|
+-------------------------------------+--------------------------+------------------------------------------------+---------------------------+
| Quantum regime                      | :math:`\chi \sim 1`      | :math:`\varepsilon_\gamma \sim \varepsilon_\pm`| Monte-Carlo               |
|                                     |                          |                                                |                           |
+-------------------------------------+--------------------------+------------------------------------------------+---------------------------+

When QED effects are negligible in the so-called classical regime (:math:`\chi \sim 10^{-3}`),
radiation back reaction can be treated as a
continuous friction force acting on the particles.
Several models have been published such as the LAD, Landau-Lifshitz ([Landau1947]_),
the model of Sokolov and of Capdessus.
The ones used in :program:`Smilei` are
based on the Landau-Lifshitz model approximated at high gamma factors
(:math:`\gamma \gg 1`).

In the quantum regime, photons with energies of the order of the energies of
the emitting electron can be produced (:math:`\varepsilon_\gamma \sim \varepsilon_\pm`).
This is treated using a Monte-Carlo
description of discrete high-energy photon emission.

In the intermediate regime (:math:`\chi \sim 1`), where the energy of the emitted photons remains
small with respect to that of the emitting electrons, but for which the
stochastic nature of photon emission cannot be neglected, the electron dynamics
is described by the addition of a stochastic term derived from a Fokker-Planck
expansion ([Niel2017]_).

Tab :numref:`radiationRegimes` can be used to well configure the radiation loss
in :program:`Smilei` (see :ref:`the radiation configuration in Species <Species>`).

The next sections describe in more details the models implemented
in :program:`Smilei` to deal with the different regimes of emission.

----

II. Continuous radiation models
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Approximated Landau-Lifshitz classical model
""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""

The classical radiation friction force acting on an ultra-relativistic electron
has been derived in [Landau1947]_. Equation of evolution of momentum that is implemented in
PIC codes is composed of both the classical Lorentz force :math:`F_L`
and the radiation friction term :math:`F_{rad}` so that:

.. math::
  :label: momentumEq

  \frac{d\mathbf{p}}{dt} = \mathbf{F}_L + \mathbf{F}_{rad}

with

.. math::
  :label: LLFrictionForce

  \mathbf{F}_{rad} = -\frac{2}{3} e \tau_e \gamma \left( \frac{d\mathbf{E}}{dt} + \mathbf{u} \times \frac{\mathbf{B}}{dt} \right) \\
  + \frac{2}{3} \frac{e}{E_{cr}} \left[ \left( \mathbf{u} \cdot \mathbf{E} \right) \mathbf{E} - \mathbf{B} \times \left( \mathbf{E} + \mathbf{u} \times \mathbf{B} \right) \right] \\
  - \frac{2}{3}\frac{e}{E_{cr}} \gamma^2 \left[ \left( \mathbf{E} + \mathbf{u} \times \mathbf{B} \right)^2 - \left( \mathbf{u} \cdot \mathbf{E}\right)^2 \right] \mathbf{u}

where :math:`\mathbf{u} = \mathbf{p} / (\gamma m c)` is the velocity,
:math:`\mathbf{p}` the momentum,
:math:`\tau_e = r_e / c = e^2 / 4 \pi \varepsilon_0 m c^3`
the time for light to travel across the classical electron radius
and :math:`E_{cr} = E_s / \alpha`
is the critical field and :math:`\alpha = e^2 / \hbar c 4 \pi  \varepsilon_0`
the fine structure constant.

For an ultra-relativistic electron, :math:`\gamma \gg 1`, some terms in
Eq. :eq:`LLFrictionForce` not explicited here can be neglected so that the
friction force reduces to a single term:

.. math::
  :label: LLFrictionForceApprox

  \mathbf{F}_{rad} = - P_{cl} \mathbf{u} / \left( \mathbf{u} c^2 \right)

where :math:`P_{cl} = \frac{2}{3} \frac{\alpha^2 mc^2}{\tau_e} \chi^2`.

The corresponding emitted power distribution as a function of the photon
frequency :math:`\omega` reads

.. math::
  :label: ClasRadPower

  \frac{dP}{d\omega} = \frac{9 \sqrt{3}}{8 \pi} \frac{P_{cl}}{ \omega_c}
  \frac{\omega}{\omega_c} \int_{\omega/\omega_c}^{+\infty}{dy K_{5/3}(y)}

with :math:`K_\nu(z)` the modified Bessel function of the second kind,
:math:`\omega_c = 3 \gamma \alpha \chi / (2 \tau_e)` the critical frequency for
synchrotron emission.
This classical approach requires the emitted photon energy
:math:`\varepsilon_\gamma = \hbar\omega` to be much smaller than that of
the emitting particle. This translates to :math:`\chi \ll 1` as given in the
introduction. Otherwise, the radiated power is know to strongly overestimate
the physical radiated energy when :math:`\chi` approaches 0.1.

Corrected classical model
""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""

In the quantum emission regime, under the conditions of slowly varying arbitrary
field compared to the formation time of the radiated photon (i) and
undercritical in respect to the Schwinger field (ii), the Lorentz invariant production
rate of high-energy photons via the multiphoton inverse Compton scattering
can be written as:

.. math::
  :label: PhotonProdRate

  \frac{d^2N}{dt d\chi_\gamma} = \frac{1}{\pi \sqrt{3}} \frac{\alpha^2}{\tau_e \chi_\pm}
  \left[ \int_\nu^{+\infty}{K_{5/3(y)}dy} + \frac{2 \chi_\gamma \nu}{2} K_{2/3}(\nu) \right]

Conditions (i) is fulfilled when :math:`a_0 = e \| A^{\mu} \| / mc^2 \gg 1`, :math:`A^{\mu}`
being the four-potential laser amplitude.
conditions (ii) corresponds to :math:`\mathbf{B}^2 - \mathbf{E}^2 \ll E_s^2`
and  :math:`\mathbf{B}\cdot \mathbf{E} \ll 1`.

From Eq. :eq:`PhotonProdRate` can be deduced the emitted power distribution in
term of the photon normalized energy. After integration, one obtains the
expression of the radiated power in the quantum regime:

.. math::
  :label: quantumRadPower

  P_{rad} = P_{cl} g(\chi_{\pm})

with

.. math::
  :label: g

  g \left( \chi_{\pm} \right) = \frac{9 \sqrt{3} }{8 \pi} \int_0^{+\infty}{d\nu
  \left[  \frac{2\nu^2 }{\left( 2 + 3 \nu \chi_\pm \right) ^2}K_{5/3}(\nu) +
  \frac{4 \nu \left( 3 \nu \chi_\pm\right)^2 }{\left( 2 + 3 \nu \chi_\pm \right)^4}K_{2/3}(\nu) \right]}

The quantum instantaneous radiated power is nothing else than the classical one
multiplied by a correction function called :math:`g \left( \chi_{\pm} \right)`.

We can simply use Eq. :eq:`LLFrictionForceApprox` with this correction close to
1 when :math:`\chi_{\pm} \ll 1` and rapidly dropping otherwise.
Thanks to this correction, the radiated energy is correct but this model does
not take into account the stochastic effects induced when the photon energy is
closed to the emitting electron. This is the subject of the next sections.

III. Stochastic schemes
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Fokker-Planck stochastic model
""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""

The Fokker-Planck approach is an extension of the corrected Landau-Lifshitz
model with an operator that takes into account diffusive stochastic effects
([Niel2017]_):

.. math::
  :label: NielStochasticForce

  F_{rad} dt = \left[ -P_{cl} g \left( \chi \right) dt + mc^2
  \sqrt{R\left( \chi, \gamma \right)} dW \right]
  \mathbf{u} / \left( \mathbf{u} c^2 \right)

where :math:`dW` is a Wiener process of variance :math:`dt`.

.. math::
  :label: NielR

    R\left( \chi, \gamma \right) = \frac{2}{3} \frac{\alpha^2}{\tau_e} \gamma
    h \left( \chi \right)

.. math::
  :label: Nielh

    h \left( \chi \right) = \frac{9 \sqrt{3}}{4 \pi} \int_0^{+\infty}{d\nu
    \left[ \frac{2\chi^3 \nu^3}{\left( 2 + 3\nu\chi \right)^3} K_{5/3}(\nu)
    + \frac{54 \chi^5 \nu^4}{\left( 2 + 3 \nu \chi \right)^5} K_{2/3}(\nu) \right]}

Monte-Carlo quantum model
""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""

The Monte-Carlo treatment of the emission is a more complex process than
the previous ones that can be divided into several steps ([Duclous2011]_,
[Lobet2013]_, [Lobet2015]_):

1. The particle (electron or positron) is first assigned a final optical depth
:math:`\tau_f` sampled from :math:`\tau_f = -\log{\xi}` where :math:`\xi` is a
random number between 0 and 1. Emission occurs when
this final optical depth is reached. A incremental optical depth :math:`\tau`
is therefore secondly assigned initially set to 0.

2. The optical depth :math:`\tau` then evolves according to the field and particle
energy variations following this integral:

.. math::
  :label: MCDtauDt

    \frac{d\tau}{dt} = \int_0^{\chi_\pm}{ \frac{d^2N}{d\chi dt}  d\chi }

that is nothing else than the production rate of photons
(integration of Eq. :eq:`PhotonProdRate`).

3. The emitted photon quantum parameter :math:`\chi_\gamma` is computed by
inverting the cumulative distribution function:

.. math::
  :label: CumulativeDistr

    P(\chi_\pm,\chi_\gamma) = \frac{\displaystyle{\int_0^{\chi_\gamma}{F(\chi_\pm, \chi)
    d\chi}}}{\displaystyle{\int_0^{\chi_\pm}{F(\chi_\pm, \chi) d\chi}}}

:math:`F` is the so-called synchrotron emissivity function so that

.. math::
  :label: MCF

    \frac{d^2 N}{dt d\chi} = \frac{2}{3} \frac{\alpha^2}{\tau_e} F (\chi, \chi_\gamma)

Inversion of  :math:`P(\chi_\pm,\chi_\gamma)=\xi'` is done after drawing
a second random number
:math:`\xi' \in \left[ 0,1\right]` to find :math:`\chi_\gamma`.

4. The energy of the emitted photon is then computed:
:math:`\varepsilon_\gamma = mc^2 \gamma_\gamma =
mc^2 \gamma_\pm \chi_\gamma / \chi_\pm`.

5. The particle momentum is then updated using momentum conservation
considering forward emission (valid when :math:`\gamma_\pm \gg 1`).

.. math::
  :label: momentumUpdate

    F_{rad} = - \frac{\varepsilon_\gamma}{c} \frac{\mathbf{p_\pm}}{\| \mathbf{p_\pm} \|}

The radiated force is just the recoil induced by the photon emission.
Radiation loss is therefore a discrete process.
Note that momentum conservation does not exactly conserve energy.
It can be shown that the error :math:`\epsilon` tends to 0 when the particle
energy tends to infinity ([Lobet2015]_) and that the error is low when
:math:`\varepsilon_\pm \gg 1` and :math:`\varepsilon_\gamma \ll \varepsilon_\pm`.
Between emission events, the electron dynamics is still governed by the
Lorentz force.

V. Implementation
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Classes for the radiations are located in the directory `src/Radiation`.
In :program:`Smilei`, the radiative process is not incorporated in the pusher.
The process is done separately using a `factory` as for the pusher and ionization.
This decision has been taken in order to:

* preserve the vector performance of the pusher when using non-vectorizable
  radiation model such as the Monte-Carlo process.
* be consistent with the current implementation
* easily be able to use any pusher (without making the code more complex)

Description of the files:

* Class `RadiationTable`: useful tools, parameters and the tables.
* Class `Radiation`: the generic one from which will inherit specific
  classes for each model.
* Class `RadiationFactory`: manage the choice of the correct radiation model
  depending of the species.
* Class `RadiationLandauLifshitz`: classical Landau-Lifshitz radiation process.
* Class `RadiationCorrLandauLifshitz`: corrected Landau-Lifshitz radiation process.
* Class `RadiationNiel`: stochastic diffusive model of [Niel2017]_.
* Class `RadiationMonteCarlo`: Monte-Carlo model.

As explained in the following, many functions have been tabulated because of
the cost of their computation for each particle. This table can be generated by
:program:`Smilei` at the initialization.
The parameters such as the ranges and the discretization can be
given in the namelist :ref:`Radiations <Radiations>`.
Once generated, the table can be written on the disk and reloaded for a next run.
The location of the table can be given in :ref:`Radiations <Radiations>`.
Small tables coded in hdf5 are provided in the repository in the folder
databases with the name: `radiation_tables.h5`.

Landau-Lifshitz based models
""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""

The classical Landau-Lifshitz model approximated for high-:math:`\gamma`
given by Eq. :eq:`LLFrictionForceApprox`
has been implemented in :program:`Smilei`
using a simple explicit scheme.
The model is accessible in the species configuration under the name
`Landau-Lifshitz`.

For the corrected version, we use a fit of the function
:math:`g(\chi)` given by Eq. :eq:`quantumCorrFit`.

.. math::
  :label: quantumCorrFit

  g \left( \chi_{\pm} \right) = \left[ 1 + 4.8 \left( 1 + \chi_{\pm} \right)
  \log \left( 1 + 1.7 \chi_{\pm} \right) + 2.44 \chi_{\pm}^2 \right]^{-2/3}

This fit enables to keep the vectorization of the particle loop.
The corrected model is accessible in the species configuration under the name
`corrected-Landau-Lifshitz`

Fokker-Planck stochastic model
""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""

Eq. :eq:`NielStochasticForce` is implemented in :program:`Smilei` using
a simple explicit scheme.

Eq. :eq:`Nielh` is tabulated for performance issue.
A polynomial fit of this integral can also be obtained in log-log
or log10-log10 domain. However, high accuracy requires high-order polynomials.
(order 20 for an accuracy around :math:`10^{-10}` for instance)

This table can be generated by :program:`Smilei` at the initialization.
The parameters such as the :math:`\chi` range and the discretization can be
given in the namelist :ref:`Radiations <Radiations>`.

The stochastic diffusive model is accessible in the species configuration
under the name `Niel`.

Monte-Carlo quantum model
""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""

The computation of Eq. :eq:`MCDtauDt` would be too expensive for every single
particles. Instead, the integral of the function :math:`F` is tabulated.
This table is referred to as `integfochi` in the code.

Similarly, Eq. :eq:`CumulativeDistr` is tabulated.
This table is referred to as `xip` in the code.
The only difference is that we choose a minimum photon quantum parameter
:math:`\chi_{\gamma,\min}` for the integration so that:

.. math::
  :label: chiMin

    \frac{\displaystyle{\int_{0}^{\chi_{\gamma,\min}}{F(\chi_\pm, \chi)
    d\chi}}}{\displaystyle{\int_0^{\chi_\pm}{F(\chi_\pm, \chi) d\chi}}} < \epsilon

This enables to find a lower bound to the :math:`\chi_\gamma` range
(discretization in the log domain) so that the
remaining part is negligible in term of radiated energy.
The parameter :math:`\epsilon` is called `xip_threshold` in the
:ref:`Radiations <Radiations>` namelist.

The tables can be generated by :program:`Smilei` at the initialization.
The parameters such as the :math:`\chi` range and the discretization can be
given in the namelist :ref:`Radiations <Radiations>`.

The Monte-Carlo model is accessible in the species configuration
under the name `Monte-Carlo`.

V. Benchmarks
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Counter-propagating Plane Wave 1D
""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""

This benchmark is referred to as `tst1d_9_rad_counter_prop.py` in the benchmark
folder. In this benchmark, a GeV electron bunch is initialized near the right
domain boundary and propagates to the left boundary from which a plane
wave is injected. The laser has an amplitude of :math:`a_0 = 270`
corresponding to an intensity of :math:`10^{23}\ \mathrm{Wcm^{-2}}` at
:math:`\lambda = 1\ \mathrm{\mu m}`. The maximal quantum parameter :math:`\chi`
value reached during the simulation is around 0.5.

.. _rad_counter_prop_scalar:

.. figure:: _static/rad_counter_prop_scalar.png
  :width: 15cm

  Comparison of the model energy scalar diagnostics. The kinetic, radiated and total
  energy are respectively plotted with solid, dashed and dotted lines for
  the Monte-Carlo (**MC**, blue), Niel (**Niel**, orange),
  corrected Landau-Lifshitz (**CLL**, green) and the Landau-Lifshitz models
  (**LL**, red).

Evolution of the kinetic, radiated and total energy is shown in
:numref:`rad_counter_prop_scalar` for all models.
The Monte-Carlo, the Niel and the corrected Landau-Lifshitz models exhibit close
results in term of total radiated and kinetic energy evolution with a final
radiation rate of 80% the initial kinetic energy. The relative error on the
total energy is small of the order of :math:`3\times10^{-3}`.
As expected, the Landau-Lifshitz (in red) overestimates the radiated energy
because the interaction happens mainly in the quantum regime.

.. _rad_counter_prop_track:

.. figure:: _static/rad_counter_prop_track.png
  :width: 18cm

  Comparison of the evolution of the normalized kinetic energy
  :math:`\gamma - 1` for some selected electrons between the radiative models
  (**CLL** for corrected Landau-Lifshitz and **LL** for Landau-Lifshitz).

:numref:`rad_counter_prop_track` shows the evolution of the normalized
kinetic energy for some selected electrons. The Monte-Carlo and the Niel models
reproduce the stochastic nature of the trajectories in comparison with the
continuous approach (corrected Landau-Lifshitz and Landau-Lifshitz).
In the latter one, every particles initially located at the same position will
follow the same trajectories.
The stochastic nature of the emission for high :math:`\chi` values can
have consequences in term of final spatial and energy distributions.
Not shown here, the Niel stochastic model do not reproduce correctly the
moment of order 3 as explained in [Niel2017]_.

Synchrotron 2D
""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""

A bunch of electron evolves in a constant magnetic field orthogonal to their
initial propagation direction.
This corresponds to two different scripts in the benchmark folder:

* `tst2d_8_synchrotron_chi1.py`: This script tests and compares the corrected
  Landau-Lifshitz and the Monte-Carlo model for an initial :math:`\chi = 1`.
* `tst2d_9_synchrotron_chi0.1.py`: This script tests and compares the corrected
  Landau-Lifshitz and the Niel model for an initial :math:`\chi = 0.1`.

In this section, we focus on the case with initial quantum parameter
:math:`\chi = 0.1`.
The magnetic field amplitude is :math:`B = 90 m \omega_r / e`.
Initial electron Lorentz factor is around :math:`\gamma =450 mc`.

Time evolution of the kinetic energy, the radiated energy and the total energy
is shown in Fig. :numref:`synchrotron_scalar`.

.. _synchrotron_scalar:

.. figure:: _static/synchrotron_scalar.png
  :width: 15cm

  Comparison of the energy scalar diagnostics. The kinetic, radiated and total
  energy are respectively plotted with solid, dashed and dotted lines for
  the Monte-Carlo (**MC**, blue), Niel (**Niel**, orange),
  corrected Landau-Lifshitz (**CLL**, green).

.. _synchrotron_x_y_gamma:

.. figure:: _static/synchrotron_x_y_gamma.png
  :width: 18cm

  Average normalized kinetic energy at simulation time :math:`25 \omega_r^{-1}`
  for the simulations with the Monte-Carlo, the Niel
  and the corrected Landau-Lifshitz (**CLL**) models.

.. _synchrotron_t_gamma_ne:

.. figure:: _static/synchrotron_t_gamma_ne.png
  :width: 18cm

  Time evolution of the electron energy distribution for the Monte-Carlo, the Niel
  and the corrected Landau-Lifshitz (**CLL**) models.

VI. Performances
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

The cost of the different models is summarized in table
:numref:`radiationTimes`.
Reported times are for the field projection, the particle pusher and
the radiation losses together. Percentages correspond to the overhead induced by
the radiation module in comparison to the standard PIC pusher.
We use different short keywords for the radiative models:

* **None**: there is no radiation loss.
* **LL**: the classical Landau-Lifshitz model approximated for large :math:`\gamma`.
* **CLL**: the Landau-Lifshitz model with the quantum correction.
* **Niel**: the stochastic model of Niel `et al.`.
* **MC**: the Monte-Carlo radiative model.

All the presented numbers are not generalizable and are only indicated to give
an idea of the model costs. the creation of macro-photons is not enabled for
the Monte-Carlo radiation process.

.. _radiationTimes:

+-------------------------------------+------------+----------+--------------+----------+--------+
| Radiation model:                    | None       | LL       | CLL          | Niel     | MC     |
+=====================================+============+==========+==============+==========+========+
| Counter-propagating Plane Wave 1D   | 0.25s      | 0.3s     | 0.36s        | 0.54s    | 0.84s  |
+-------------------------------------+------------+----------+--------------+----------+--------+
| Synchrotron 2D                      | 3.9s       | 4.2s     | 4.8s         | 9s       | 5.6s   |
| :math:`\chi=0.5`,  :math:`B=100`    |            | - 10%    | - 30%        | - 140%   | - 50%  |
+-------------------------------------+------------+----------+--------------+----------+--------+
| Interaction with a carbon thin foil | 6.5s       | 6.9s     | 7.2s         | 7.4s     | 7.2s   |
| 2D                                  |            |          |              |          |        |
+-------------------------------------+------------+----------+--------------+----------+--------+

Descriptions of the cases:

* **Counter-propagating Plane Wave 1D**: Collision between an electron bunch
  and a counter-propagating plane wave.
  The case is run on a single node of Jureca with 2 MPI ranks and 12 OpenMP
  threads per rank.

* **Synchrotron 2D**: The domain is fulfilled with electrons having the same
  initial momentum so that initially :math:`\chi=0.5` with the constant magnetic
  field :math:`B_z=100`. The domain has a dimension of 496x496 cells with
  16 particles per cell and 8x8 patches.
  A 4th order B-spline shape factor is used for the projection.
  The case is run on a single node of Jureca with 2 MPI ranks and 12 OpenMP
  threads per rank.

* **Thin foil 2D**:
  This case simulates the interaction of a fully-ionized carbon thin foil
  with an extremely intense plane wave in 2D.
  The thin foil in located at 4 :math:`\mu\mathrm{m}` of the left border `xmin`.
  It starts with a linear preplasma of 1 :math:`\mu\mathrm{m}` followed with
  a uniform section of 3 :math:`\mu\mathrm{m}` of density 492 :math:`n_c`.
  The target is irradiated by a Gaussian plane wave of peak intensity
  :math:`a_0 = 270` corresponding to :math:`10^{23}\ \mathrm{Wcm^{-2}}` and FWHM 50 fs.
  The domain has a discretization of 64 cells per :math:`\mu\mathrm{m}` in
  the two directions x and y with 64 particles per cell.
  There is two species: electrons and carbon ions.
  Only electrons can radiate.
  The case is run on 16 nodes of Poincare with 2 MPI ranks and 8 OpenMP
  threads per rank.

For the moment, only LL and CLL can be vectorized efficiently
as for the pushers.
As a consequence, code performance is likely to be more impacted running on
SIMD architecture with large vector registers such as Intel Xeon Phi.
