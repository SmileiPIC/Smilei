Synchrotron-like radiation loss
-------------------------------

High-energy particles traveling in a strong electromagnetic field radiates
their energy. Depending on the field strength and the particle energy, radiation
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
continuous friction force acting on the particles based on the Landau-Lifshitz
model approximated at high gamma factors.

In the quantum regime, photons with energies of the order of the energies of
the emitting electron can be produced. This is treated using a Monte-Carlo
description of discrete high-energy photon emission.

----

Physical model
^^^^^^^^^^^^^^

The quantum Lorentz invariant parameter :math:`\chi_\pm` is an indicator of
the radiation  regime,

.. math::
  :label: particleQuantumParameter

  \chi_{\pm} = \frac{\gamma_{\pm}}{E_s} \left| \left( \beta \cdot \mathbf{E}
  \right)^2 - \left( \mathbf{E} + \mathbf{v} \times \mathbf{B} \right)^2
  \right|^{1/2}

where :math:`E_s = 4 \pi \alpha \varepsilon_0 m^2 c^4 / e^3` is the Schwinger field.

Continuous classical schemes
^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Approximated Landau-Lifshitz classical model
""""""""""""""""""""""""""""""""""""""""""""

Corrected classical model
"""""""""""""""""""""""""

Stochastic semi-classical or quantum schemes
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Fokker-Planck model
""""""""""""""""""

Monte-Carlo quantum model
"""""""""""""""""""""""""


Benchmarks
^^^^^^^^^^

Performances
^^^^^^^^^^^^
