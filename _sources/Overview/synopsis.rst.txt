Synopsis
--------

Smilei is a collaborative project providing physicists with an open-source,
user-friendly, high-performance and multi-purpose
**electromagnetic Particle-In-Cell (PIC) code for plasma simulation**.

The code is developed in C++ based on an object-oriented architecture.
To face the diverse needs of the Smilei community, it offers modularity:

* various geometries (Cartesian 1D, 2D, 3D or cylindrical with decomposition into azimuthal modes),
* arbitrary laser or plasma profiles (any Python function), 
* various Maxwell solvers, particle pushers, interpolators, projectors
* an envelope solver, including in the cylindrical geometry
* advanced boundary conditions (e.g. Perfectly-Matched Layers)
* etc.

The **user-friendly interface** consists in input files written in the Python language,
and a whole set of run-time diagnostics (outputs in HDF5) and user-friendly (Python)
:doc:`post-processing </Use/post-processing>` tools complement the code.

Co-developed by HPC specialists and physicists, Smilei is **designed for high performances**
on massively-parallel super-computers. It benefits from a state-of-the-art hybrid
MPI/OpenMP parallelization, dynamic load balancing and SIMD vectorization.
It has been successfully tested on various architectures, among which the most recent
Intel Cascadelake (CSL) & Fujitsu A64FX (ARM).

Recently, GPU acceleration has been implemented in SMILEI and allows offloading on Nvidia or AMD GPUs, such as V100, A100 or MI250.
As of yet, not all features are supported.  

In Smilei, Maxwell’s equations are solved using a Yee mesh, where the
electric and magnetic fields are centered following the finite-difference time-domain (FDTD)
method or related methods. A pseudo-spectral analytical time domain method is
also available as an experimental feature.
Charge deposition follows a charge-conservation scheme.

As a multi-purpose code, Smilei is applied to a **wide range of physics-related studies**:
from relativistic laser-plasma interaction to astrophysics. Smilei thus benefits from
various additional physics modules, among which :ref:`field ionization <field_ionization>`,
:doc:`binary collisions and impact ionization </Understand/collisions>`. QED processes, such as
high-energy photon emission and its :doc:`back-reaction </Understand/radiation_loss>`
on the electron dynamics, as well as
:doc:`pair production </Understand/multiphoton_Breit_Wheeler>` through the Breit-Wheeler
process, are also included.

An initial detailed account (as of 2018) of Smilei's capabilities is given in
`this article <https://doi.org/10.1016/j.cpc.2017.09.024>`_.
For publications on more advanced features, please refer to the :doc:`material` section of this documentation. 
