Releases
--------

Major releases are available here as well as on the
`GitHub page <https://github.com/SmileiPIC/Smilei>`_.
We greatly appreciate external users trying this code and giving feedback.

You can :doc:`contact us <partners>` to become a developer of the official releases:
the `developer's Gitlab repository <https://llrgit.in2p3.fr/smilei/smilei>`_ is used
for development, and you will need to ask for a password first,
but the process is fast and straightforward.

----

Upcoming changes
^^^^^^^^^^^^^^^^

* Monte Carlo QED photon emission
* Moving window along y and z
* Tracking particles in *Screens*

----

.. _latestVersion:

Current release 3.2
^^^^^^^^^^^^^^^^^^^

**Download**: `Smilei v3.2 <_downloads/smilei-v3.2.tar.gz>`_

* New pushers (Vay's and Higuera-Cary's)
* *Numpy* used for filtering track particles
* Fourth order in 3D
* Add some missing 3D features: external fields management, boundary conditions and non-neutral plasma initialization
* OpenMP support in moving window
* Tracked particles post-processing improved for large files
* Bugfixes: energy computation in 3D or with moving window, random number seed

----

Release 3.1
^^^^^^^^^^^

**Download**: `Smilei v3.1 <_downloads/smilei-v3.1.tar.gz>`_

* *Screen* diagnostics
* Exporting 3D diagnostics to VTK for reading in ParaView or VisIt
* Partial support of the `OpenPMD <https://www.openpmd.org>`_ standard
* Improvements: moving window (OpenMP), 3D projection
* Bugfixes: tracked particles, walls, collisional ionization, etc.

Notes:

* Outputs of Fields and Tracks are incompatible with 3.0
* The input "output_dir" is not supported anymore

----

Release 3.0
^^^^^^^^^^^

**Download**: `Smilei v3.0 <_downloads/smilei-v3.0.tar.gz>`_

* **3D geometry**
* Field and scalar diagnostics improved for more flexibility and memory saving
* Faster initialization (including Maxwell-Jüttner sampling)
* Post-processing handles restarts
* Bugfixes in checkpoints, timers, memory profile

----

Release 2.3
^^^^^^^^^^^

**Download**: `Smilei v2.3 <_downloads/smilei-v2.3.tar.gz>`_

* Post-processing scripts have been turned into a *python* module
* Many bugfixes, such as addressing diagnostics efficiency


----

Release 2.0
^^^^^^^^^^^

**Download**: `Smilei v2.2 <_downloads/smilei-v2.2.tar.gz>`_

* **state-of-the-art dynamic load balancing**
* full *python* namelist, allowing for complex, user-friendly input
* external fields and antennas
* binary Coulomb collisions
* new diagnostics
* *python* scripts for post-processing

----

Release 1.0
^^^^^^^^^^^

**Download**: `Smilei v1.0 <_downloads/smilei-v1.0.tar.gz>`_

* 1D & 2D cartesian geometries
* Moving window
* Hybrid MPI-OpenMP parallelization
* Field ionization
* Some python diagnostics

