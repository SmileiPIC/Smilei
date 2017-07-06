Synchrotron-like radiation loss
-------------------------------

I. Introduction
^^^^^^^^^^^^^^^

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

Different approaches have been implemented in :program:`Smilei`.
When QED effects are negligible in the so-called classical regime,
radiation back reaction can be treated as a
continuous friction force acting on the particles.
Several models have been published. The ones used in :program:`Smilei` are
based on the Landau-Lifshitz model approximated at high gamma factors.

In the quantum regime, photons with energies of the order of the energies of
the emitting electron can be produced. This is treated using a Monte-Carlo
description of discrete high-energy photon emission.

In the intermediate regime, where the energy of the emitted photons remains
small with respect to that of the emitting electrons, but for which the
stochastic nature of photon emission cannot be neglected, the electron dynamics
is described by the addition of a stochastic term derived from a Fokker-Planck
expansion.

----

II. Physical models
^^^^^^^^^^^^^^^^^^^

The quantum Lorentz invariant parameter :math:`\chi_\pm` is an indicator of
the radiation  regime,

.. math::
  :label: particleQuantumParameter

  \chi_{\pm} = \frac{\gamma_{\pm}}{E_s} \left| \left( \beta \cdot \mathbf{E}
  \right)^2 - \left( \mathbf{E} + \mathbf{v} \times \mathbf{B} \right)^2
  \right|^{1/2}

where :math:`E_s = m^2 c^3 / \hbar e \simeq 1.3 10^{18}\ \mathrm{V/m}` is
the Schwinger field, :math:`e` is the electron charge,
:math:`m` is the electron mass, :math:`c` the speed of light in vacuum,
:math:`\hbar` the reduced Planck constant.

:math:`\alpha = e^2 / \hbar c 4 \pi  \varepsilon_0` the fine structure constant

III. Implementations
^^^^^^^^^^^^^^^^^^^^

Continuous radiation models
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Approximated Landau-Lifshitz classical model
""""""""""""""""""""""""""""""""""""""""""""

Corrected classical model
"""""""""""""""""""""""""

Stochastic schemes
^^^^^^^^^^^^^^^^^^

Fokker-Planck stochastic model
""""""""""""""""""""""""""""""

Monte-Carlo quantum model
"""""""""""""""""""""""""


IV. Benchmarks
^^^^^^^^^^^^^^

V. Performances
^^^^^^^^^^^^^^^

The cost of the different models is summarized in table
:numref:`radiationTimes`.
Reported times are for the field projection, the particle pusher and
the radiation losses together.

.. _radiationTimes:

+-------------------------------------+------------+----------+--------------+----------+--------+
| Radiation model:                    | None       | LL       | Corrected LL | Niel     | MC     |
+=====================================+============+==========+==============+==========+========+
| Synchrotron 2D                      | 3.9s       | 4.2s     | 4.8s         | 9s       | 5.6s   |
| :math:`\chi=0.5`,  :math:`B=100`    |            |          |              |          |        |
+-------------------------------------+------------+----------+--------------+----------+--------+

Descriptions of the cases:

* **Synchrotron 2D**: the domain is fulfilled with electrons having the same
  initial momentum so that initially :math:`\chi=0.5` with the constant magnetic field
  :math:`B_z=100`. The domain has a dimension of 128x128 cells with
  16 particles per cell and 8x8 patches.
  A 4th order B-spline shape factor is used for the projection.
  The case is run on a single node of Jureca with 2 MPI ranks and 12 OpenMP
  threads per rank.
