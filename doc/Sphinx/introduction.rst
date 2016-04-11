Introduction
------------

To face the diverse needs of the teams involved in its development, :program:`Smilei`
is developed in C++, based on an object-oriented architecture.
Its modularity allows to run simulations in various dimensions, geometries, to chose
between various Maxwell solvers, particle pushers, interpolators, projectors, etc.
Today, the one-dimensional in space, three-dimensional in velocity (1D3V)
and 2D3V versions of the code have been developed and benchmarked.

Maxwell's equations are solved on the so-called Yee-mesh with centered electric
and magnetic fields using the finite-difference time-domain (FDTD) method
or related methods [Nuter2014]_\ . Moreover, charge deposition is computed
following a charge-conservation scheme [Esirkepov2001]_\ . 

:doc:`Binary collisions <collisions>` have been implemented and
Monte-Carlo routines are currently under development to account for
(i) high-energy (gamma) photon emission and its back-reaction on the electron 
dynamics, as well as (ii) electron-positron pair creation. These developments are
undertaken in collaboration with the team that has introduced similar routines
in the PIC code :program:`Calder` (see [Lobet2013]_). Such routines will be of
particular importance for the modelling of strongly relativistic astrophysical scenarii.

On the performance side, :program:`Smilei` benefits from a state-of-the-art
hybrid **MPI/OpenMP parallelization**, an original particle sorting algorithm,
and innovative dynamic load balancing.
:program:`Smilei` is therefore designed to run on massively parallel machines,
and its flexibility should allow one to benefit from the newest and futures HPC architectures.

----

.. _latestVersion:

Current release 2.0
^^^^^^^^^^^^^^^^^^^

**Download**: `Smilei v2.0 (pre-release) <_downloads/smilei-v2.0.tar.gz>`_

A second release is **in preparation** featuring great improvements:

* **state-of-the-art** dynamic load balancing
* full *python* namelist, allowing for complex, user-friendly input
* 3D cartesian geometry
* external fields and antennas
* binary collisions
* new diagnostics
* *python* scripts for post-processing
* various improvements on stability and error management

We greatly appreciate external users trying this code and giving feedback.
At the moment, developer rights are still restricted,
but you can :doc:`contact us <partners>` to become a developer.


----

Release 1.0
^^^^^^^^^^^

**Download**: `Smilei v1.0 <_downloads/smilei-v1.0.tar.gz>`_

After nearly two years of development, :program:`Smilei` v1.0 offers some nice features:

* 1D & 2D cartesian geometries
* moving-window
* hybrid MPI-OpenMP parallelization
* field ionization
* some python diagnostics
* open-source

It is protected by a "licence CeCILL", the french equivalent to the Gnu GPL license
for "logiciels libres".

.. warning::
  This version does not have associated documentation.
  Please refer to the examples.



----

Example : Electron Acceleration
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Below, an example of electron acceleration by laser wakefield.
The figure represents the evolution of the electronic density in time. 
A hotspot of electron is created behind the bubble.

.. raw:: html

    <video controls="controls">
    <source src="_static/Rho_electron1long.ogg" type="video/ogg" />
    </video>

----

Scalability
^^^^^^^^^^^

.. rst-class:: inprogress
  
  In progress ...

.. rubric :: 1. OpenMP: Electron Acceleration

The hotspot of electrons produces an important imbalance between the
compute load of the different MPI processes involved in the simulation.

OpenMP permits to smooth this phenomenon by balancing macro-particles between threads.

.. image:: _static/perfsOMP.png
    :width: 500px


.. rubric :: 2. MPI: SBS Amplification

In the completely opposite context of a very homogeneous plasma, we oberve during a
"Grand challenge" on `Occigen <https://www.cines.fr/calcul/materiels/occigen>`_,
a good scaling at very large scale.

.. image:: _static/SMILEI_Scaling.png
    :width: 500px

----

References
^^^^^^^^^^

.. [Nuter2014] Nuter *et al.*, Eur. J. Phys. D **68**, Issue 6 (2014)

.. [Esirkepov2001] Esirkepov, Comp. Phys. Comm. **135**, 144 (2001)

.. [Lobet2013] Lobet *et al.*, arXiv:1311.1107 (2013), Plasma Phys. Control. Fusion
  




