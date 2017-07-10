Synchrotron-like radiation loss
-------------------------------

I. Introduction
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

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
|                                     |                          |                                                |                           |
+-------------------------------------+--------------------------+------------------------------------------------+---------------------------+
| Semi-classical radiation emission   | :math:`\chi \in 10^{-3}  |                                                | Corrected Landau-Lifshitz |
|                                     | - 10^{-2}`               |                                                |                           |
+-------------------------------------+--------------------------+------------------------------------------------+---------------------------+
| Weak quantum regime                 | :math:`\chi \sim 10^{-1}`|                                                | Stochastic model of       |
|                                     |                          |                                                | Niel `et al` / Monte-Carlo|
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
expansion.

Tab :numref:`radiationRegimes` can be used to well configure the radiation loss
in :program:`Smilei` (see :ref:`the radiation configuration in Species <Species>`).

The next sections describe in more details the models implemented
in :program:`Smilei` to deal with the different regimes of emission.

----

II. Continuous radiation models
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Approximated Landau-Lifshitz classical model
""""""""""""""""""""""""""""""""""""""""""""

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

where :math:`\mathbf{u} = \mathbf{p} / (\gamma m c)` is the normalized momentum,
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

Eq. :eq:`LLFrictionForceApprox` has been implemented in :program:`Smilei`
under the radiation model name `Landau-Lifshitz`.

Corrected classical model
"""""""""""""""""""""""""

III. Stochastic schemes
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Fokker-Planck stochastic model
""""""""""""""""""""""""""""""

Monte-Carlo quantum model
"""""""""""""""""""""""""


IV. Benchmarks
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

V. Performances
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

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
| Interaction with a carbon thin foil | 6.5s       | 6.9s     | 7.2s         | 7.8s     | 7.2s   |
| 2D                                  |            |          |              |          |        |
+-------------------------------------+------------+----------+--------------+----------+--------+

Descriptions of the cases:

* **Counter-propagating Plane Wave 1D**: Collision between an electron bunch
  and a counter-propagating plane wave.
  The case is run on a single node of Poincare with 2 MPI ranks and 8 OpenMP
  threads per rank.

* **Synchrotron 2D**: The domain is fulfilled with electrons having the same
  initial momentum so that initially :math:`\chi=0.5` with the constant magnetic
  field :math:`B_z=100`. The domain has a dimension of 128x128 cells with
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
