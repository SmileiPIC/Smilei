Releases
--------

Major releases are available here as well as on the
`GitHub page <https://github.com/SmileiPIC/Smilei>`_.
We greatly appreciate external users trying this code and giving feedback.
You can submit *issues* when experiencing difficulties,
or *pull requests* for your changes to become part of the official releases.

Note that most of the development of the code is currently hosted in
a `different repository <https://llrgit.in2p3.fr/smilei/smilei>`_
reserved for the :doc:`partners`. It is regularly synchronized with
the GitHub page.

.. warning::

  In v3.3, :doc:`significant changes<syntax_changes>` have been made to the input syntax.

----

Upcoming changes
^^^^^^^^^^^^^^^^

* Vectorization
* Interface with the PICSAR library
* Faster collisions
* Probes can output species' current and density

----

.. _latestVersion:

Release 3.5
^^^^^^^^^^^^^^^^^^^^^

**Download**: `Smilei v3.5 <_downloads/smilei-v3.5.tar.gz>`_

* :doc:`Laser defined in tilted plane<laser_offset>`
* Bugfixes: Field diagnostic subgrid, Scalar diagnostic PoyInst, MPI tags for large number of patches

----

Release 3.4.1
^^^^^^^^^^^^^^^^^^^^^

**Download**: `Smilei v3.4.1 <_downloads/smilei-v3.4.1.tar.gz>`_

* Ionization considering a user-defined rate

----

Release 3.4
^^^^^^^^^^^

**Download**: `Smilei v3.4 <_downloads/smilei-v3.4.tar.gz>`_

* Compatibility with Python 3
* New 'Performances' diagnostic
* Tracked particles may output the fields at their location
* 'subgrid' option in Fields diagnostics
* Printout of the expected disk usage
* Laser propagation pre-processing
* More flexible domain decomposition
* Relativistic initialization
* Particles injection using Numpy arrays
* Possibility to use user-defined ionization rates
* Bugfixes: circular polarization, collisional ionization

----

Release 3.3
^^^^^^^^^^^

**Download**: `Smilei v3.3 <_downloads/smilei-v3.3.tar.gz>`_

* **Major** :doc:`syntax changes<syntax_changes>` in the namelist
* QED radiation reaction
* Monte-Carlo QED photon emission
* *Test mode* to quickly check the namelist consistency
* *ParticleBinning* and *Screen* diagnostics accept a python function as their ``deposited_quantity`` and ``axis``.
* Bugfixes: 4th order, field ionization

----

Release 3.2
^^^^^^^^^^^

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

Release 2.2
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
