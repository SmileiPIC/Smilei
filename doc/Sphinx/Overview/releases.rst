Releases
--------


Get Smilei
^^^^^^^^^^^^^^^^

**Clone the latest version of Smilei from** `GitHub <https://github.com/SmileiPIC/Smilei>`_:

.. code ::

  git clone https://github.com/SmileiPIC/Smilei.git

*Learn about Git* `here <https://git-scm.com/doc>`_.


You can find older, unsupported versions here <https://github.com/SmileiPIC/Smilei/releases>

----

.. _latestVersion:

Changes made in the repository (not released)
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

* Much faster ``DiagFields`` (speedup ~ x3)
* Collisions: new parameter ``time_frozen``
* Performances diagnostic: new parameter ``cumulative``
* Laser Envelope: multi-level tunnel ionization creates multiple electrons, improving the sampling
* For developers: new table management for Monte-Carlo physical processes (transparent to users)
* Bugfixes: 

  * Poisson Solver correction was not properly accounted for with SDMD.
  * Bug correction using Monte-Carlo radiation and multiphoton Breit-Wheeler processes with checkpoints
  * C++11 compilation issue
  * Reading particle weights and momenta from hdf5 file
  * solved segfault with Multiphoton Breit-Wheeler process in `AMcylindrical` geometry


----

Projects
^^^^^^^^^^^^^^^^

* Already available, but experimental:

  * :doc:`/Understand/SDMD`
  * Particle merging
  * Nuclear reactions
  * Interface with the PICSAR library for AM spectral solver

* In preparation:

  * Perfectly-matched layers for the envelope model
  * More spectral solvers
  * GPU support


----

Release 4.7
^^^^^^^^^^^^^^^^^^^^^

* **Perfectly Matched Layers** boundary conditions for EM fields (+2D Cartesian benchmark).
* Improved performance for ARM-based processors including the Fujitsu A64FX
* Improved performance for GNU, LLVM, arm-clang and Fujitsu compilers on all types of architectures
* Lasers can be injected from all boundaries
* Flag ``ponderomotive_dynamics`` removed from ``Species`` block. All ``Species`` interact with ``LaserEnvelope`` if present 
* Option to create neutrons for D-D fusion
* Collisions can be done less often
* Lasers can be injected from all boundaries
* New 4th-order non-standard FDTD solver ``M4``
* Timestep dependent field interpolation scheme
* ``LaserOffset``:

  * may be re-used from a previous simulation
  * available from ``ymin``, ``ymax``, ``zmin`` and ``zmax``
  * has new arguments ``fft_time_window`` and ``fft_time_step``

* Diagnostics:

  * Probes can include components of the Poynting vector ``PoyX``, ``PoyY``, ``PoyZ``
  * Probes can be time-integrated
  * ``ParticleBinning`` diagnostics may accept ``"auto"`` as axis limits
  * Particle IDs may be modified in the ``DiagTrackParticles.filter`` (8 available bits)
  * Screens may have a ``cylinder`` shape
  * Scalar diagnostics for AM geometry now available
  * happi ``ParticleBinning`` now uses the keyword ``average`` instead of ``sum``

* Bugfixes:

  * Poynting scalars behaviour with several patches, or with checkpoints
  * Densities too low are put to 0 to avoid underflow
  * Prescribed fields in 2D
  * ``ellipticity = -1.`` was doing ``+1.``
  * Setting ``units`` in happi's ``TrackParticles`` was wrong (for plotting only)
  * Current communication correction for FIR filters
  * Fix for particle merging segmentation fault in spherical and Cartesian modes
  * Tracked particles with the vectorized mode
  * ``momentum_initialization`` from a file did not take the proper file

----

Release 4.6
^^^^^^^^^^^^^^^^^^^^^

* :doc:`/Understand/SDMD`
* New 4th-order non-standard FDTD solver ``Bouchard`` for 2D and 3D geometries
* New method for current filtering with a user-provided FIR kernel for 1D, 2D and 3D geometries
* Diagnostics may now have a ``name`` (useful during post-processing)
* Laser Envelope:

  * linear and circular polarization
  * ionization model
  * normalized laser frequency can be different from 1

* Particles can be imported from a file
* Some :doc:`/Use/profiles` can be imported from a file
* Coulomb logarithm may be multiplied by a constant factor
* Happi:

  * handles fonts
  * time slider available with multiple plotting
  * ``vsym`` option for symmetric graph
  * ``getXmoved`` now accounts for requested units
  * Tracked particles can be selected before sorting

* Bugfixes:

  * Fix in the vectorized projection at order 4
  * Photons could not be read from numpy array
  * DiagFields with ``time_average`` did not work for densities
  * Prescribed fields caused unstable real fields
  * Initialisation from numpy or hdf5 caused wrong weights in AM geometry
  * Better positionning of collisionally-ionised electrons
  * Fix segfault from thermalizing boundary
  * Running a simulation displayed the wrong version v4.4

----

Release 4.5
^^^^^^^^^^^^^^^^^^^^^

* Changes:

  * Current filtering with adjustable number of passes per dimension
  * Improved axial boundary conditions for ``AMcylindrical`` geometry
  * Units in ``RadiationSpectrum`` diagnostic are more consistent with that
    of ``ParticleBinning``
  * Ionisation current at fourth order of interpolation
  * Correction for :doc:`/Understand/collisions` as suggested in [Higginson2020]_

* Bugfixes:

  * ``PrescribedField`` was sometimes not applied by some OpenMP threads
  * Scalar ``Ukin_bnd`` was sometimes wrong with load balancing
  * Scalar ``Urad`` was sometimes wrong with moving window
  * On some systems, particles IDs were incorrect with ionization


----

Release 4.4
^^^^^^^^^^^^^^^^^^^^^

* Changed radiation tables: see :doc:`the doc </Understand/radiation_loss>`.

  * :red:`Old tables are not valid anymore, input files must be updated.`
  * Default tables are now embebded in the code
  * Possibility to read external generated by an :doc:`external tool </Use/tables>` (more efficient and stable)

* New ``RadiationSpectrum`` diagnostics available (see :doc:`the doc </Understand/radiation_loss>`)
* ``AMcylindrical``: sorting, documentation, subgrid in DiagFields,
  species-related currents and density in probes (not per mode anymore)
* LaserOffset is not recomputed after restart
* Prescribed fields that only contribute to pushing particles
* Laser Envelope: added envelope equation solver with reduced numerical dispersion
* Bugfixes:

  * Weight-initialization bug in AM geometry when a species was initialized
    on top of a regularly-initialized species
  * LaserOffset was off sideways and temporally by a couple of cells
  * Do not project twice a frozen species
  * Probes for species faulty when 4th order of interpolation
  * Checkpoints ``restart_number=0`` was not used
  * Checkpointing with ``dump_minutes`` could be out of sync between MPI process
  * Prevent deadlock when restart files are corrupted
  * Checkpoints ``file_grouping`` had typo with python3
  * Scalar ``Ukin`` for ions was incorrect, thus ``Ubal`` was also wrong
  * happi had incorrect unit conversion with a sum of two fields
  * fix error occurring when envelope Probes on axis are used in AM geometry


----

Release 4.3
^^^^^^^^^^^^^^^^^^^^^


* ``AMcylindrical`` : envelope, ionization, additional diagnotics,
  number of ppc per direction, binomial current filter, poisson solver,
  non-separable laser initialization per mode, improved diag field nomenclature
* Particle injector
* More control over the moving window movement
* More control over the regular position initialization in Cartesian geometries
* Bugfixes:

  * ionization of frozen species
  * particle binning was not following the moving window
  * gaussian profile with order 0 was incorrect
  * tracked particles post-processing was incorrect above 20M particles
  * better management of particle binning in collisions
  * Intel 19 optimizations


----

Release 4.2
^^^^^^^^^^^^^^^^^^^^^

* ``AMcylindrical`` geometry with azimuthal Fourier decomposition (beta version)
* Different convention for circular polarization amplitude
* 1D and 2D laser envelope model
* Compatibility between various ionization and QED models
* Bugfixes:

  * Binomial filter in Cartesian 3D parallel implementation
  * Various crashes linked to vectorization
  * ``LaserGaussian2D`` when focused far from boundary
  * Laser :py:data:`a0` normalization to :py:data:`omega`
  * Frozen particles are now properly ionized
  * Position initialization over another species with moving window
  * Tracked particles output was missing the mass factor for momenta
  * Breit-Wheeler pair production with fine grain sorted particles


----

Release 4.1
^^^^^^^^^^^^^^^^^^^^^

* Probe diagnostics of currents and density per species
* Field diagnostics with more than 2^32 points
* Bugfixes:

  * collisions (badly affected by vectorization)
  * adaptive vectorization with dynamic load balancing
  * memory leak in the laser envelope model

* Disable usage of ``-ipo`` to compile on supercomputers
  despite of saving time simulation

  * it needs too many resources (time and memory) to link
  * it is recommended to do some tests on a new supercomputer
    without and then to re-establish it

.. warning::

  Since version 4.1, the :ref:`definition of macro-particle weights<Weights>`
  has changed to ensure they do not depend on the cell volume. This impacts
  only the users working directly with values of weights. Other simulation
  results should be unchanged.


----

Release 4.0
^^^^^^^^^^^^^^^^^^^^^

* :ref:`vectorization`
* :ref:`laser_envelope`
* MPI option ``MPI_THREAD_MULTIPLE`` is now optional (but recommended)
* Faster collisions
* Bugfixes: handling ``sum`` for happi's ``ParticleBinning``

----

Release 3.5
^^^^^^^^^^^^^^^^^^^^^

* :doc:`Laser defined in tilted plane</Use/laser_offset>`
* Bugfixes: Field diagnostic subgrid, Scalar diagnostic PoyInst,
  MPI tags for large number of patches

----

Release 3.4.1
^^^^^^^^^^^^^^^^^^^^^

* Ionization considering a user-defined rate

----

Release 3.4
^^^^^^^^^^^

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

* **Major** :doc:`syntax changes</syntax_changes>` in the namelist
* QED radiation reaction
* Monte-Carlo QED photon emission
* *Test mode* to quickly check the namelist consistency
* *ParticleBinning* and *Screen* diagnostics accept a python function as their
  ``deposited_quantity`` and ``axis``.
* Bugfixes: 4th order, field ionization

----

Release 3.2
^^^^^^^^^^^

* New pushers (Vay's and Higuera-Cary's)
* *Numpy* used for filtering track particles
* Fourth order in 3D
* Add some missing 3D features: external fields management, boundary conditions
  and non-neutral plasma initialization
* OpenMP support in moving window
* Tracked particles post-processing improved for large files
* Bugfixes: energy computation in 3D or with moving window, random number seed

----

Release 3.1
^^^^^^^^^^^

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

* **3D geometry**
* Field and scalar diagnostics improved for more flexibility and memory saving
* Faster initialization (including Maxwell-Jüttner sampling)
* Post-processing handles restarts
* Bugfixes in checkpoints, timers, memory profile

----

Release 2.3
^^^^^^^^^^^

* Post-processing scripts have been turned into a *python* module
* Many bugfixes, such as addressing diagnostics efficiency


----

Release 2.2
^^^^^^^^^^^

* **state-of-the-art dynamic load balancing**
* full *python* namelist, allowing for complex, user-friendly input
* external fields and antennas
* binary Coulomb collisions
* new diagnostics
* *python* scripts for post-processing

----

Release 1.0
^^^^^^^^^^^

* 1D & 2D cartesian geometries
* Moving window
* Hybrid MPI-OpenMP parallelization
* Field ionization
* Some python diagnostics
