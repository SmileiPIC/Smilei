Synopsis
--------

To face the diverse needs of the teams involved in its development, :program:`Smilei`
is developed in C++, based on an object-oriented architecture.
Its modularity allows to run simulations in various dimensions, geometries, to chose
between various Maxwell solvers, particle pushers, interpolators, projectors, etc.
:program:`Smilei` works now in 3D cartesian geometry.

Maxwell's equations are solved on the so-called Yee-mesh with centered electric
and magnetic fields using the finite-difference time-domain (FDTD) method
or related methods. Moreover, charge deposition is computed
following a charge-conservation scheme.

:doc:`Binary collisions <collisions>` have been implemented and
Monte-Carlo routines are currently under development to account for
high-energy (gamma) photon emission and its back-reaction on the electron 
dynamics, as well as electron-positron pair creation. These developments are
undertaken in collaboration with the team that introduced
`similar routines <http://iopscience.iop.org/article/10.1088/1742-6596/688/1/012058>`_
in the PIC code :program:`Calder`. Such routines will be of
particular importance for the modelling of strongly relativistic astrophysical scenarii.

On the performance side, :program:`Smilei` benefits from a state-of-the-art
hybrid **MPI/OpenMP parallelization**, an original particle sorting algorithm,
and innovative dynamic load balancing.
:program:`Smilei` is therefore designed to run on massively parallel machines,
and its flexibility should allow one to benefit from the newest and futures HPC architectures.

