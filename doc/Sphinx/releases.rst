Releases
--------

This page lists the major changes, but it is recommended to
get the latest version of Smilei on `GitHub <https://github.com/SmileiPIC/Smilei>`_.

----

Upcoming changes
^^^^^^^^^^^^^^^^

* Interface with the PICSAR library (currently experimental)
* Particle merging (beta version)


----

.. _latestVersion:

Latest version
^^^^^^^^^^^^^^^^^^^^^

The latest version tarball can be donwloaded here:

**Download**: `Smilei latest <_downloads/Smilei.tar.gz>`_


----

Release 4.2
^^^^^^^^^^^^^^^^^^^^^

**Download**: `Smilei v4.2 <_downloads/smilei-v4.2.tar.gz>`_


* Different convention for circular polarization amplitude
* Binomial filter in Cartesian 3D bug fix in parallel implementation
* 1D and 2D laser envelope model
* Cylindrical geometry with azimuthal Fourier decomposition (beta version)
* Compatibility between various ionization and QED models
* Bugfixes:
   * Various crashes linked to vectorization
   * `LaserGaussian2D` when focused far from boundary
   * Laser :py:data:`a0` normalization to :py:data:`omega`
   * Frozen particles are now properly ionized
   * Position initialization over another species with moving window
   * Tracked particles output was missing the mass factor for momenta
   * Breit-Wheeler pair production with fine grain sorted particles


----

Release 4.1
^^^^^^^^^^^^^^^^^^^^^

**Download**: `Smilei v4.1 <_downloads/smilei-v4.1.tar.gz>`_

* Probe diagnostics of currents and density per species
* Field diagnostics with more than 2^32 points
* Bugfixes:

 * collisions (badly affected by vectorization)
 * adaptive vectorization with dynamic load balancing
 * memory leak in the laser envelope model
 
* Disable usage of `-ipo` to compile on supercomputers despite of saving time simulation

 * it needs too many resources (time and memory) to link
 * it is recommended to do some tests on a new supercomputer without and then to re-establish it

.. warning::

  Since version 4.1, the :ref:`definition of macro-particle weights<Weights>`
  has changed to ensure they do not depend on the cell volume. This impacts
  only the users working directly with values of weights. Other simulation
  results should be unchanged.


----

Release 4.0
^^^^^^^^^^^^^^^^^^^^^

**Download**: `Smilei v4.0 <_downloads/smilei-v4.0.tar.gz>`_

* :ref:`vectorization`
* :ref:`laser_envelope`
* MPI option `MPI_THREAD_MULTIPLE` is now optional (but recommended)
* Faster collisions
* Bugfixes: handling `sum` for happi's `ParticleBinning`

----

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
