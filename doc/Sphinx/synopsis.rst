Synopsis
--------

Smilei is a collaborative project providing the scientific community with an open-source,
user-friendly, high-performance and multi-purpose Particle-In-Cell (PIC) code
for plasma simulation.

To face the diverse needs of the Smilei community, the code is developed in C++,
based on an object-oriented architecture. Its modularity and user-friendly Python
interface allows to run simulations in various dimensions (Cartesian 1D, 2D, 3D),
to chose between various Maxwell solvers, particle pushers, interpolators, projectors, etc.
A whole set of run-time diagnostics (based on HDF5) and user-friendly (Python)
:doc:`post-processing <post-processing>` tools complements the code.

Co-developed by HPC specialists and physicists, Smilei is designed for high performances
on massively-parallel super-computers. It benefits from a state-of-the-art hybrid
MPI/OpenMP parallelization, dynamic load balancing and SIMD vectorization.
It has been successfully tested on various architectures among which the most recent
Intel Skylake (SKL) & Knigths Landing (KNL).

In Smilei, Maxwell’s equations are solved using a Yee mesh with centered
electric and magnetic fields following the finite-difference time-domain (FDTD)
method or related methods. A pseudo-spectral analytical time domain method is
also available as an experimental feature.
Charge deposition follows a charge-conservation scheme.

As a multi-purpose code, Smilei is applied to a wide range of physics-related studies:
from relativistic laser-plasma interaction to astrophysics. Smilei thus benefits from
various additional physics modules among which :ref:`field ionization <field_ionization>`,
:doc:`binary collisions and impact ionization <collisions>`. QED processes, such as
high-energy photon emission and its :doc:`back-reaction <radiation_loss>`
on the electron dynamics, as well as
:doc:`pair production <multiphoton_Breit_Wheeler>` through the Breit-Wheeler
process, are also included.

A detailed account of Smilei's capabilities is given in
`this article <https://doi.org/10.1016/j.cpc.2017.09.024>`_.
