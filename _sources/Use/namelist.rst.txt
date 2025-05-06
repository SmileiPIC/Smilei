.. |exp| raw:: html

   <span class="exp-label">experimental</span>

Write a namelist
----------------

Before you run :program:`Smilei`, you need a *namelist* (an input file). The namelist
is written in the *python* language. It is thus recommended to know the basics of *python*.

We suggest you copy one existing namelist from the folder *benchmarks*.
All namelists have the extension ``.py``.


----

General rules
^^^^^^^^^^^^^

* :program:`Smilei` requires a few *blocks* to be defined, such as::

    Main(
        # ...
        timestep = 0.01,         # defines the timestep value
        grid_length = [10., 20.], # defines the 2D box dimensions
        # ...
    )

  Outside blocks, you can calculate anything you require.
  Inside a block, you must only define variables for :program:`Smilei`.

* The *python* syntax requires special indentation of each line.
  You begin with no indentation, but you have to **add four spaces at the
  beginning of lines inside a group**, and so on.
  For instance::

    if a == 0:
        timestep = 0.1
        if b == 1:
            timestep = 0.2
    else:
        timestep = 0.3

* You will need to use `lists <https://docs.python.org/2/tutorial/introduction.html#lists>`_,
  which are series of things in *python*,
  defined between brackets ``[]`` and separated by commas.
  For example, ``mean_velocity = [0., 1.1, 3.]``.

* You are free to import any installed *python* package into the namelist.
  For instance, you may obtain :math:`\pi` using ``from math import pi``.

* All quantities are normalized to arbitrary values: see :doc:`/Understand/units`.

----

Python workflow
^^^^^^^^^^^^^^^

*Python* is started at the beginning of the simulation (one *python* interpreter
for each MPI process). The following steps are executed:

#. A few variables from :program:`Smilei` are passed to *python* so that they are
   available to the user:

   * The rank of the current MPI process as :py:data:`smilei_mpi_rank`.
   * The total number of MPI processes as :py:data:`smilei_mpi_size`.
   * The number of OpenMP threads per MPI :py:data:`smilei_omp_threads`.
   * The total number of cores :py:data:`smilei_total_cores`.

#. The namelist(s) is executed.

#. *Python* runs :py:data:`preprocess()` if the user has defined it.
   This is a good place to calculate things that are not needed for
   post-processing with :program:`happi`.

#. The simulation is initialized (including field and particle arrays).

#. *Python* runs :py:data:`cleanup()` if the user has defined it.
   This is a good place to delete unused heavy variables.

#. *Python* checks whether the *python* interpreter is needed during the simulation
   (e.g. the user has defined a temporal :doc:`profile <profiles>` which requires *python*
   to calculate it every timestep). Otherwise, *python* is stopped.

All these instructions are summarized in a file ``smilei.py``,
so that the user can directly run ``python -i smilei.py`` for post-processing purposes.

----

Main variables
^^^^^^^^^^^^^^

The block ``Main`` is **mandatory** and has the following syntax::

  Main(
      geometry = "1Dcartesian",
      interpolation_order = 2,
      interpolator = "momentum-conserving",
      grid_length  = [16. ],
      cell_length = [0.01],
      simulation_time    = 15.,
      timestep    = 0.005,
      number_of_patches = [64],
      cluster_width = 5,
      maxwell_solver = 'Yee',
      EM_boundary_conditions = [
          ["silver-muller", "silver-muller"],
  #        ["silver-muller", "silver-muller"],
  #        ["silver-muller", "silver-muller"],
      ],
      time_fields_frozen = 0.,
      reference_angular_frequency_SI = 0.,
      print_every = 100,
      random_seed = 0,
  )

.. py:data:: geometry

  The geometry of the simulation:

  * ``"1Dcartesian"``
  * ``"2Dcartesian"``
  * ``"3Dcartesian"``
  * ``"AMcylindrical"``: cylindrical geometry with :doc:`/Understand/azimuthal_modes_decomposition`.

  In the following documentation, all references to dimensions or coordinates
  depend on the ``geometry``.
  1D, 2D and 3D stand for 1-dimensional, 2-dimensional and 3-dimensional cartesian
  geometries, respectively. All coordinates are ordered as :math:`(x)`, :math:`(x,y)` or :math:`(x,y,z)`.
  In the ``"AMcylindrical"`` case, all grid coordinates are 2-dimensional
  :math:`(x,r)`, while particle coordinates (in :ref:`Species`)
  are expressed in the 3-dimensional Cartesian frame :math:`(x,y,z)`.

  .. warning::

    The ``"AMcylindrical"`` geometry has some restrictions.
    Boundary conditions must be set to ``"remove"`` for particles,
    ``"silver-muller"`` for longitudinal EM boundaries and
    ``"buneman"`` for transverse EM boundaries.
    You can alternatively use ``"PML"`` for any EM boundary.
    Collisions and
    order-4 interpolation are not supported yet.

.. py:data:: interpolation_order

  :default: ``2``

  Interpolation order, defines particle shape function:

  * ``1``  : 2 points stencil in r with Ruyten correction, 3 points stencil in x. Supported only in AM geometry.
  * ``2``  : 3 points stencil, supported in all configurations.
  * ``4``  : 5 points stencil, not supported in vectorized 2D geometry.

  The Ruyten correction is the scheme described bu equation 4.2 in `this paper <https://www.sciencedirect.com/science/article/abs/pii/S0021999183710703>`_ .
  It allows for a more accurate description on axis at the cost of a higher statistic noise so it often requires the use of more macro-particles.

.. py:data:: interpolator

  :default: ``"momentum-conserving"``

  * ``"momentum-conserving"``
  * ``"wt"``

  The interpolation scheme to be used in the simulation.
  ``"wt"`` is for the timestep dependent field interpolation scheme described in
  `this paper <https://doi.org/10.1016/j.jcp.2020.109388>`_ .

.. py:data:: grid_length
             number_of_cells

  A list of numbers: size of the simulation box for each dimension of the simulation.

  * Either ``grid_length``, the simulation length in each direction in units of :math:`L_r`,
  * or ``number_of_cells``, the number of cells in each direction.

  .. note::
    
    In ``AMcylindrical`` geometry, the grid represents 2-dimensional fields.
    The second dimension is the **radius** of the cylinder.

.. py:data:: cell_length

  A list of floats: sizes of one cell in each direction in units of :math:`L_r`.


.. py:data:: simulation_time
             number_of_timesteps

  Duration of the simulation.

  * Either ``simulation_time``, the simulation duration in units of :math:`T_r`,
  * or ``number_of_timesteps``, the total number of timesteps.


.. py:data:: timestep
             timestep_over_CFL

  Duration of one timestep.

  * Either ``timestep``, in units of :math:`T_r`,
  * or ``timestep_over_CFL``, in units of the *Courant–Friedrichs–Lewy* (CFL) time.

.. py:data:: gpu_computing

   :default: ``False``
   
   Activates GPU acceleration if set to True

.. py:data:: number_of_patches

  A list of integers: the number of patches in each direction.
  Each integer must be a power of 2, and the total number of patches must be
  greater or equal than the number of MPI processes.
  It is also strongly advised to have more patches than the total number of openMP threads.
  See :doc:`/Understand/parallelization`.On the other hand, in case of GPU-acceleration it is recommended to use one patch per MPI-rank
  (with one MPI-rank per GPU)


.. py:data:: patch_arrangement

  :default: ``"hilbertian"``

  Determines the ordering of patches and the way they are separated into the
  various MPI processes. Options are:

  * ``"hilbertian"``: following the Hilbert curve (see :ref:`this explanation<LoadBalancingExplanation>`).
  * ``"linearized_XY"`` in 2D or ``"linearized_XYZ"`` in 3D: following the
    row-major (C-style) ordering.
  * ``"linearized_YX"`` in 2D or ``"linearized_ZYX"`` in 3D: following the
    column-major (fortran-style) ordering. This prevents the usage of
    :ref:`Fields diagnostics<DiagFields>` (see :doc:`/Understand/parallelization`).

.. py:data:: cluster_width

  :default: set to minimize the memory footprint of the particles pusher, especially interpolation and projection processes

  For advanced users. Integer specifying the cluster width along X direction in number of cells.
  The "cluster" is a sub-patch structure in which particles are sorted for cache improvement.
  ``cluster_width`` must divide the number of cells in one patch (in dimension X).
  The finest sorting is achieved with ``cluster_width=1`` and no sorting with ``cluster_width`` equal to the full size of a patch along dimension X.
  The cluster size in dimension Y and Z is always the full extent of the patch.

.. py:data:: maxwell_solver

  :default: 'Yee'

  The solver for Maxwell's equations.
  Only ``"Yee"`` and ``"M4"`` are available for all geometries at the moment.
  ``"Cowan"``, ``"Grassi"``, ``"Lehe"`` and ``"Bouchard"`` are available for ``2DCartesian``.
  ``"Lehe"`` and ``"Bouchard"`` are available for ``3DCartesian``.
  ``"Lehe"`` and ``"Terzani"`` are available for ``AMcylindrical``.
  The M4 solver is described in `this paper <https://doi.org/10.1016/j.jcp.2020.109388>`_.
  The Lehe solver is described in `this paper <https://journals.aps.org/prab/abstract/10.1103/PhysRevSTAB.16.021301>`_.
  The Bouchard solver is described in `this thesis p. 109 <https://tel.archives-ouvertes.fr/tel-02967252>`_.
  The Terzani solver is described in `this paper <https://doi.org/10.1016/j.cpc.2019.04.007>`_.

.. py:data:: solve_poisson

   :default: True

   Decides if Poisson correction must be applied or not initially.

.. py:data:: poisson_max_iteration

  :default: 50000

  Maximum number of iteration for the Poisson solver.

.. py:data:: poisson_max_error

  :default: 1e-14

  Maximum error for the Poisson solver.

.. py:data:: solve_relativistic_poisson

   :default: False

   Decides if relativistic Poisson problem must be solved for at least one species.
   See :doc:`/Understand/relativistic_fields_initialization` for more details.

.. py:data:: relativistic_poisson_max_iteration

  :default: 50000

  Maximum number of iteration for the Poisson solver.

.. py:data:: relativistic_poisson_max_error

  :default: 1e-22

  Maximum error for the Poisson solver.

.. py:data:: EM_boundary_conditions

  :type: list of lists of strings
  :default: ``[["periodic"]]``

  The boundary conditions for the electromagnetic fields. Each boundary may have one of
  the following conditions: ``"periodic"``, ``"silver-muller"``, ``"reflective"``, ``"ramp??"`` or ``"PML"``.

  | **Syntax 1:** ``[[bc_all]]``, identical for all boundaries.
  | **Syntax 2:** ``[[bc_X], [bc_Y], ...]``, different depending on x, y or z.
  | **Syntax 3:** ``[[bc_Xmin, bc_Xmax], ...]``,  different on each boundary.

  * ``"silver-muller"`` is an open boundary condition.
    The incident wave vector :math:`k_{inc}` on each face is defined by
    ``"EM_boundary_conditions_k"``.
    When using ``"silver-muller"`` as an injecting boundary,
    make sure :math:`k_{inc}` is aligned with the wave you are injecting.
    When using ``"silver-muller"`` as an absorbing boundary,
    the optimal wave absorption on a given face will be along :math:`k_{abs}`
    the specular reflection of :math:`k_{inc}` on the considered face.

  * ``"ramp??"`` is a basic, open boundary condition designed
    for the spectral solver in ``AMcylindrical`` geometry.
    The ``??`` is an integer representing a number of cells
    (smaller than the number of ghost cells).
    Over the first half, the fields remain untouched.
    Over the second half, all fields are progressively reduced down to zero.

  * ``"PML"`` stands for Perfectly Matched Layer. It is an open boundary condition.
    The number of cells in the layer must be defined by ``"number_of_pml_cells"``.
    It supports laser injection as in ``"silver-muller"``.
    If not all boundary conditions are ``PML``, make sure to set ``number_of_pml_cells=0`` on boundaries not using PML.

.. py:data:: EM_boundary_conditions_k

  :type: list of lists of floats
  :default: ``[[1.,0.],[-1.,0.],[0.,1.],[0.,-1.]]`` in 2D
  :default: ``[[1.,0.,0.],[-1.,0.,0.],[0.,1.,0.],[0.,-1.,0.],[0.,0.,1.],[0.,0.,-1.]]`` in 3D

  For ``silver-muller`` absorbing boundaries,
  the *x,y,z* coordinates of the unit wave vector ``k`` incident on each face
  (sequentially Xmin, Xmax, Ymin, Ymax, Zmin, Zmax).
  The number of coordinates is equal to the dimension of the simulation.
  The number of given vectors must be equal to 1 or to the number of faces
  which is twice the dimension of the simulation. In cylindrical geometry,
  ``k`` coordinates are given in the ``xr`` frame and only the Rmax face is affected.

  | **Syntax 1:** ``[[1,0,0]]``, identical for all boundaries.
  | **Syntax 2:** ``[[1,0,0],[-1,0,0], ...]``,  different on each boundary.

.. py:data:: number_of_pml_cells

  :type: List of lists of integers
  :default: ``[[10,10],[10,10],[10,10]]``

  Defines the number of cells in the ``"PML"`` layers using the same alternative syntaxes as ``"EM_boundary_conditions"``.

.. rst-class:: experimental

.. py:data:: pml_sigma

  :type: List of profiles
  :default: [lambda x : 20 * x**2]

  Defines the sigma profiles across the transverse dimension of the PML for each dimension of the simulation.
  It must be expressed as a list of profiles (1 per dimension).

  If a single profile is given, it will be used for all dimensions.

  For a given dimension, the same profile is applied to both sides of the domain.

  The profile is given as a single variable function defined on the interval [0,1] where 0 is the inner bound of the PML and 1 is the outer bound of the PML. 
  Please refer to :doc:`/Understand/PML` if needed in AM geometry.

.. rst-class:: experimental

.. py:data:: pml_kappa

  :type: List of profiles
  :default: [lambda x : 1 + 79 * x**4]

  Defines the kappa profiles across the transverse dimension of the PML for each dimension of the simulation.
  It must be expressed as a list of profiles (1 per dimension).

  If a single profile is given, it will be used for all dimensions.

  For a given dimension, the same profile is applied to both sides of the domain.

  The profile is given as a single variable function defined on the interval [0,1] where 0 is the inner bound of the PML and 1 is the outer bound of the PML. 
  Please refer to :doc:`/Understand/PML` if needed in AM geometry.

.. py:data:: time_fields_frozen

  :default: 0.

  Time, at the beginning of the simulation, during which fields are frozen.


.. _reference_angular_frequency_SI:

.. py:data:: reference_angular_frequency_SI

  The value of the reference angular frequency :math:`\omega_r` in SI units,
  **only needed when collisions, ionization, radiation losses
  or multiphoton Breit-Wheeler pair creation are requested**.
  This frequency is related to the normalization length according to :math:`L_r\omega_r = c`
  (see :doc:`/Understand/units`).


.. py:data:: print_every

  Number of timesteps between each info output on screen. By default, 10 outputs per
  simulation.


.. py:data:: print_expected_disk_usage

  :default: ``True``

  If ``False``, the calculation of the expected disk usage, that is usually printed in the
  standard output, is skipped. This might be useful in rare cases where this calculation
  is costly.


.. py:data:: random_seed

  :default: 0

  The value of the random seed. Each patch has its own random number generator, with a seed
  equal to ``random_seed`` + the index of the patch.

.. py:data:: number_of_AM

  :type: integer
  :default: 2

  The number of azimuthal modes used for the Fourier decomposition in ``"AMcylindrical"`` geometry.
  The modes range from mode 0 to mode ``"number_of_AM-1"``.

.. py:data:: number_of_AM_classical_Poisson_solver

  :default: 1

  The number of azimuthal modes used for the field initialization with non relativistic Poisson solver in ``"AMcylindrical"`` geometry.
  Note that this number must be lower or equal to the number of modes of the simulation.

.. py:data:: number_of_AM_relativistic_field_initialization

  :default: 1

  The number of azimuthal modes used for the relativistic field initialization in ``"AMcylindrical"`` geometry.
  Note that this number must be lower or equal to the number of modes of the simulation.

.. py:data:: use_BTIS3_interpolation

  :default: ``False``

  If ``True``, the B-translated interpolation scheme 3 (or B-TIS3) described in :doc:`/Understand/algorithms` is used.

.. py:data:: custom_oversize

   :type: integer
   :default: 2

   The number of ghost-cell for each patches. The default value is set accordingly with
   the ``interpolation_order`` value.

..
  .. py:data:: spectral_solver_order

    :type: A list of integers
    :default: ``[0,0]`` in AM geometry.

    The order of the spectral solver in each dimension. Set order to zero for infinite order.
    In AM geometry, only infinite order is supported along the radial dimension.

..
  .. py:data:: initial_rotational_cleaning

    :default: ``False``

    If ``True``, uses the picsar library to do the rotational cleaning.

    Rotational cleaning corrects field initialization in spectral space
    in order to make sure that the fields at :math:`t=0` are a valid solution
    of the Maxwell equation.
    This operation is only supported in AM geometry and with picsar
    spectral solver. It requires a FFT of the full domain on a single MPI
    process so very large simulations may face problems with this procedure.

..
  .. py:data:: cell_sorting

    :default: ``False``

    If ``True``, forces the use of cell sorting for particles. This flag is
    automatically set to true if any feature requiring cell sorting is requested
    (vectorization, collisions or
    particle merging) so it is mainly a convenience for developers.

----

Load Balancing
^^^^^^^^^^^^^^

Load balancing (explained :ref:`here <LoadBalancingExplanation>`) consists in exchanging
patches (domains of the simulation box) between MPI processes to reduce the
computational load imbalance.
The block ``LoadBalancing`` is optional. If you do not define it, load balancing will
occur every 150 iterations.

.. code-block:: python

  LoadBalancing(
      initial_balance = True,
      every = 150,
      cell_load = 1.,
      frozen_particle_load = 0.1
  )

.. py:data:: initial_balance

  :default: True

  Decides if the load must be balanced at initialization. If not, the same amount of
  patches will be attributed to each MPI rank.

.. py:data:: every

  :default: 150

  Number of timesteps between each load balancing **or** a :ref:`time selection <TimeSelections>`.
  The value ``0`` suppresses all load balancing.

.. py:data:: cell_load

  :default: 1.

  Computational load of a single grid cell considered by the dynamic load balancing algorithm.
  This load is normalized to the load of a single particle.

.. py:data:: frozen_particle_load

  :default: 0.1

  Computational load of a single frozen particle considered by the dynamic load balancing algorithm.
  This load is normalized to the load of a single particle.

----

.. rst-class:: experimental

Multiple decomposition of the domain
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

The block ``MultipleDecomposition`` is necessary for spectral solvers and optional in all other cases.
When present, it activates
the :doc:`/Understand/SDMD` (SDMD) technique
which separates the decomposition of the field grids from that of the particles.
Fields are set on large sub-domain called *regions* (1 region per MPI process) while
particles are kept as small *patches* as in the standard decomposition (many patches per MPI process).
Benefits of this option are illustrated `in this paper <https://hal.archives-ouvertes.fr/hal-02973139>`_.


.. code-block:: python

  MultipleDecomposition(
      region_ghost_cells = 2
  )

.. py:data:: region_ghost_cells

   :type: integer
   :default: 2

   The number of ghost cells for each region.
   The default value is set accordingly with the ``interpolation_order``.
   The same number of ghost cells is used in all dimensions except for spectral solver in AM geometry for which the number of radial ghost cells is always automatically set to be the same as patches.


----

.. _Vectorization:

Vectorization
^^^^^^^^^^^^^^^^^^^^^

The block ``Vectorization`` is optional.
It controls the SIMD operations that can enhance the performance of some computations.
The technique is detailed in Ref. [Beck2019]_ and summarized in :doc:`this doc </Understand/vectorization>`.
It requires :ref:`additional compilation options<vectorization_flags>` to be actived.

.. code-block:: python

  Vectorization(
      mode = "adaptive",
      reconfigure_every = 20,
      initial_mode = "on"
  )

.. py:data:: mode

  :default: ``"off"``

  * ``"off"``: non-vectorized operators are used.
    Recommended when the number of particles per cell stays below 10.
  * ``"on"``: vectorized operators are used.
    Recommended when the number of particles per cell stays above 10.
    Particles are sorted per cell.
  * ``"adaptive"``: the best operators (scalar or vectorized)
    are determined and configured dynamically and locally
    (per patch and per species). For the moment this mode is only supported in ``3Dcartesian`` geometry.
    Particles are sorted per cell.

  In the ``"adaptive"`` mode, :py:data:`cluster_width` is set to the maximum.

.. py:data:: reconfigure_every

  :default: 20

  The number of timesteps between each dynamic reconfiguration of
  the vectorized operators, when using the  ``"adaptive"`` vectorization mode.
  It may be set to a :ref:`time selection <TimeSelections>` as well.


.. py:data:: initial_mode

  :default: ``off``

  Default state when the ``"adaptive"`` mode is activated
  and no particle is present in the patch.


----

.. _movingWindow:

Moving window
^^^^^^^^^^^^^

The simulated domain can move relatively to its the initial position. The "moving window"
is (almost) periodically shifted in the ``x_max`` direction.
Each "shift" consists in removing a column of patches from the ``x_min`` border and
adding a new one after the ``x_max`` border, thus changing the physical domain that the
simulation represents but keeping the same box size. This is particularly useful to
*follow* waves or plasma moving at high speed.
The frequency of the shifts is adjusted so that the average displacement velocity
over many shifts matches the velocity given by the user.
The user may ask for a given number of additional shifts at a given time.
These additional shifts are not taken into account for the evaluation of the average
velocity of the moving window.

The block ``MovingWindow`` is optional. The window does not move it you do not define it.

.. warning::

  When the window starts moving, all laser injections via Silver-Muller boundary conditions
  are immediately stopped for physical correctness.

.. code-block:: python

  MovingWindow(
      time_start = 0.,
      velocity_x = 1.,
      number_of_additional_shifts = 0.,
      additional_shifts_time = 0.,
  )


.. py:data:: time_start

  :type: Float.
  :default: 0.

  The time at which the window starts moving.


.. py:data:: velocity_x

  :type: Float.
  :default: 0.

  The average velocity of the moving window in the ``x_max`` direction. It muste be between 0 and 1.

.. py:data:: number_of_additional_shifts

  :type: Integer.
  :default: 0.

  The number of additional shifts of the moving window.

.. py:data:: additional_shifts_time

  :type: Float.
  :default: 0.

  The time at which the additional shifts are done.


.. note::

  The :ref:`particle binning diagnostics <DiagParticleBinning>` accept an "axis" called ``moving_x``
  corresponding to the ``x`` coordinate corrected by the moving window's current movement.

----

.. _CurrentFilter:

Current filtering
^^^^^^^^^^^^^^^^^

The present version of :program:`Smilei` provides a
:ref:`multi-pass binomial filter <multipassBinomialFilter>` on the current densities,
which parameters are controlled in the following block::

  CurrentFilter(
      model = "binomial",
      passes = [0],
      kernelFIR = [0.25,0.5,0.25]
  )

.. py:data:: model

  :default: ``"binomial"``

  The model for current filtering.

  * ``"binomial"`` for a binomial filter.
  * ``"customFIR"`` for a custom FIR kernel.

.. py:data:: passes

  :type: A python list of integers.
  :default: ``[0]``

  The number of passes (at each timestep) given for each dimension.
  If the list is of length 1, the same number of passes is assumed for all dimensions.

.. py:data:: kernelFIR

  :default: ``"[0.25,0.5,0.25]"``

  The FIR kernel for the ``"customFIR"`` model. The number of coefficients
  must be less than twice the number of ghost cells
  (adjusted using :py:data:`custom_oversize`).


----

.. _FieldFilter:

Field filtering
^^^^^^^^^^^^^^^^^

The present version of :program:`Smilei` provides a method for field filtering
(at the moment, only the :ref:`Friedman electric field time-filter <EfieldFilter>` is available)
which parameters are controlled in the following block::

  FieldFilter(
      model = "Friedman",
      theta = 0.,
  )

.. py:data:: model

  :default: ``"Friedman"``

  The model for field filtering. Presently, only ``"Friedman"`` field filtering is available.

.. py:data:: theta

  :default: ``0.``

  The :math:`\theta` parameter (between 0 and 1) of Friedman's method.


----

.. _Species:

Species
^^^^^^^

Each species has to be defined in a ``Species`` block::

  Species(
      name      = "electrons1",
      position_initialization = "random",
      momentum_initialization = "maxwell-juettner",
      regular_number = [],
      particles_per_cell = 100,
      mass = 1.,
      atomic_number = None,
      #maximum_charge_state = None,
      number_density = 10.,
      # charge_density = None,
      charge = -1.,
      mean_velocity = [0.],
      #mean_velocity_AM = [0.],
      temperature = [1e-10],
      boundary_conditions = [
          ["reflective", "reflective"],
      #    ["periodic", "periodic"],
      #    ["periodic", "periodic"],
      ],
      # thermal_boundary_temperature = None,
      # thermal_boundary_velocity = None,
      time_frozen = 0.0,
      # ionization_model = "none",
      # ionization_electrons = None,
      # ionization_rate = None,
      is_test = False,
      pusher = "boris",

      # Radiation reaction, for particles only:
      radiation_model = "none",
      radiation_photon_species = "photon",
      radiation_photon_sampling = 1,
      radiation_photon_gamma_threshold = 2,
      radiation_max_emissions = 10,

      # Relativistic field initialization:
      relativistic_field_initialization = "False",

      # For photon species only:
      multiphoton_Breit_Wheeler = ["electron","positron"],
      multiphoton_Breit_Wheeler_sampling = [1,1]

      # Merging
      merging_method = "vranic_spherical",
      merge_every = 5,
      merge_min_particles_per_cell = 16,
      merge_max_packet_size = 4,
      merge_min_packet_size = 4,
      merge_momentum_cell_size = [16,16,16],
  )

.. py:data:: name

  The name you want to give to this species.
  It should be more than one character and can not start with ``"m_"``.

.. py:data:: position_initialization

   The method for initialization of particle positions. Options are:

   * ``"regular"`` for regularly spaced. See :py:data:`regular_number`.
   * ``"random"`` for randomly distributed.
   * ``"centered"`` for centered in each cell (not supported in ``AMcylindrical`` geometry.
   * The :py:data:`name` of another species from which the positions are copied.
     The *source* species must have positions initialized using one of the three
     other options above, and must be defined before this species.
   * A *numpy* array or an *HDF5* file defining all the positions of the particles.
     In this case you must also provide the weight of each particle (see :ref:`Weights`).
     See :doc:`particle_initialization`.


.. py:data:: regular_number

   :type: A list of as many integers as the simulation dimension

   When ``position_initialization = "regular"``, this sets the number of evenly-spaced
   particles per cell in each direction: ``[Nx, Ny, Nz]`` in cartesian geometries and
   ``[Nx, Nr, Ntheta]`` in ``AMcylindrical`` in which case we recommend
   ``Ntheta`` :math:`\geq 4\times (` ``number_of_AM`` :math:`-1)`.
   If unset, ``particles_per_cell`` must be a power of the simulation dimension,
   for instance, a power of 2 in ``2Dcartesian``.

.. py:data:: momentum_initialization

  The method for initialization of particle momenta. Options are:

  * ``"maxwell-juettner"`` for a relativistic maxwellian (see :doc:`how it is done<maxwell-juttner>`)
  * ``"rectangular"`` for a rectangular distribution
  * ``"cold"`` for zero temperature
  * A *numpy* array or an *HDF5* file defining all the momenta of the particles.
    See :doc:`particle_initialization`.

  The first 2 distributions depend on the parameter :py:data:`temperature` explained below.

.. py:data:: particles_per_cell

  :type: float or :doc:`profile <profiles>`

  The number of particles per cell.


.. py:data:: mass

  The mass of particles, in units of the electron mass :math:`m_e`.


.. py:data:: atomic_number

  :default: 0

  The atomic number of the particles (must be below 101).
  It is required for ionization and nuclear reactions.
  It has an effect on collisions by accounting for the atomic screening
  (if not defined, or set to 0 for ions, screening is discarded as if
  the ion was fully ionized).

.. py:data:: maximum_charge_state

  :default: 0

  The maximum charge state of a species for which the ionization model is ``"from_rate"``.

.. py:data:: number_density
             charge_density

  :type: float or :doc:`profile <profiles>`

  The absolute value of the charge density or number density (choose one only)
  of the particle distribution, in units of the reference density :math:`N_r` (see :doc:`/Understand/units`).


.. py:data:: charge

  :type: float or :doc:`profile <profiles>`

  The particle charge, in units of the elementary charge :math:`e`.


.. py:data:: mean_velocity

  :type: a list of 3 floats or :doc:`profiles <profiles>`

  The initial drift velocity of the particles, in units of the speed of light :math:`c`, in the `x`, `y` and `z` directions.

  **WARNING**: For massless particles, this is actually the momentum in units of :math:`m_e c`.
  
.. py:data:: mean_velocity_AM

  :type: a list of 3 floats or :doc:`profiles <profiles>`

  The initial drift velocity of the particles, in units of the speed of light :math:`c`, in the longitudinal, radial and azimuthal directions.
  This entry is available only in ``AMcylindrical`` velocity and cannot be used if also ``mean_velocity`` is used in the same ``Species``: only one of the two can be chosen.

  **WARNING**: For massless particles, this is actually the momentum in units of :math:`m_e c`.

  **WARNING**: The initial cylindrical drift velocity is applied to each particle, thus it can be computationally demanding. 

.. py:data:: temperature

  :type: a list of 3 floats or :doc:`profiles <profiles>`
  :default: ``1e-10``

  The initial temperature of the particles, in units of :math:`m_ec^2`.


.. py:data:: boundary_conditions

  :type: a list of lists of strings
  :default: ``[["periodic"]]``

  The boundary conditions for the particles of this species.
  Each boundary may have one of the following conditions:
  ``"periodic"``, ``"reflective"``, ``"remove"`` (particles are deleted),
  ``"stop"`` (particle momenta are set to 0), and ``"thermalize"``.
  For photon species (``mass=0``), the last two options are not available.

  | **Syntax 1:** ``[[bc_all]]``, identical for all boundaries.
  | **Syntax 2:** ``[[bc_X], [bc_Y], ...]``, different depending on x, y or z.
  | **Syntax 3:** ``[[bc_Xmin, bc_Xmax], ...]``,  different on each boundary.

.. py:data:: thermal_boundary_temperature

  :default: None

  A list of floats representing the temperature of the thermal boundaries (those set to
  ``"thermalize"`` in  :py:data:`boundary_conditions`) for each spatial coordinate.
  Currently, only the first coordinate (x) is taken into account.

.. py:data:: thermal_boundary_velocity

  :default: []

  A list of floats representing the components of the particles' drift velocity after
  encountering the thermal boundaries (those set to ``"thermalize"`` in :py:data:`boundary_conditions`).

.. py:data:: time_frozen

  :default: 0.

  The time during which the particles are "frozen", in units of :math:`T_r`.
  Frozen particles do not move and therefore do not deposit any current density either.
  Nonetheless, they deposit a charge density.
  They are computationally much cheaper than non-frozen particles and oblivious to any EM-fields
  in the simulation. Note that frozen particles can be ionized (this is computationally much cheaper
  if ion motion is not relevant).

.. py:data:: ionization_model

  :default: ``"none"``

  The model for :ref:`field ionization <field_ionization>`:

  * ``"tunnel"`` for tunnel ionization using :ref:`PPT-ADK <ppt_adk>` (requires species with an :py:data:`atomic_number`)
  * ``"tunnel_full_PPT"`` |exp| for tunnel ionization using :ref:`PPT-ADK with account for magnetic number<ppt_adk>` (requires species with an :py:data:`atomic_number`)
  * ``"tunnel_envelope_averaged"`` for :ref:`field ionization with a laser envelope <field_ionization_envelope>`
  * ``"from_rate"``, relying on a :ref:`user-defined ionization rate <rate_ionization>` (requires species with a :py:data:`maximum_charge_state`).

.. py:data:: bsi_model

  :default: ``"none"``

  Apply the :ref:`Barrier Suppression Ionization <barrier_suppression>` correction for ionization in strong fields.
  This correction is supported only for ``ionization_model`` = ``"tunnel"`` or ``tunnel_full_PPT``.
  The available BSI models are:

  * ``"Tong_Lin"`` for :ref:`Tong and Lin <tong_lin>`'s rate.
  * ``"KAG"`` for :ref:`Kostyukov Artemenko Golovanov <KAG>`'s rate. 

.. py:data:: ionization_rate

  A python function giving the user-defined ionisation rate as a function of various particle attributes.
  To use this option, the `numpy package <http://www.numpy.org/>`_ must be available in your python installation.
  The function must have one argument, that you may call, for instance, ``particles``.
  This object has several attributes ``x``, ``y``, ``z``, ``px``, ``py``, ``pz``, ``charge``, ``weight`` and ``id``.
  Each of these attributes are provided as **numpy** arrays where each cell corresponds to one particle.

  The following example defines, for a species with maximum charge state of 2,
  an ionization rate that depends on the initial particle charge
  and linear in the x coordinate:

  .. code-block:: python

    from numpy import exp, zeros_like

    def my_rate(particles):
        rate = zeros_like(particles.x)
        charge_0 = (particles.charge==0)
        charge_1 = (particles.charge==1)
        rate[charge_0] = r0 * particles.x[charge_0]
        rate[charge_1] = r1 * particles.x[charge_1]
        return rate

    Species( ..., ionization_rate = my_rate )

.. py:data:: ionization_electrons

  The name of the electron species that :py:data:`ionization_model` uses when creating new electrons.


.. py:data:: is_test

  :default: ``False``

  Flag for test particles. If ``True``, this species will contain only test particles
  which do not participate in the charge and currents.


.. .. py:data:: c_part_max
..
..   :red:`to do`
..

.. py:data:: pusher

  :default: ``"boris"``

  Type of pusher to be used for this species. Options are:

  * ``"boris"``: The relativistic Boris pusher
  * ``"borisnr"``: The non-relativistic Boris pusher
  * ``"vay"``: The relativistic pusher of J. L. Vay
  * ``"higueracary"``: The relativistic pusher of A. V. Higuera and J. R. Cary
  * ``"norm"``:  For photon species only (rectilinear propagation)
  * ``"ponderomotive_boris"``: modified relativistic Boris pusher for species interacting with the laser envelope model. Valid only if the species has non-zero mass
  * ``"borisBTIS3"``: as ``"boris"``, but using B fields interpolated with the B-TIS3 scheme.
  * ``"ponderomotive_borisBTIS3"``: as ``"ponderomotive_boris"``, but using B fields interpolated with the B-TIS3 scheme.

  **WARNING**: ``"borisBTIS3"`` and ``"ponderomotive_borisBTIS3"`` can be used only when ``use_BTIS3_interpolation=True`` in the ``Main`` block.
  
.. py:data:: radiation_model

  :default: ``"none"``

  The **radiation reaction** model used for this species (see :doc:`/Understand/radiation_loss`).

  * ``"none"``: no radiation
  * ``"Landau-Lifshitz"`` (or ``ll``): Landau-Lifshitz model approximated for high energies
  * ``"corrected-Landau-Lifshitz"`` (or ``cll``): with quantum correction
  * ``"Niel"``: a `stochastic radiation model <https://arxiv.org/abs/1707.02618>`_ based on the work of Niel `et al.`.
  * ``"Monte-Carlo"`` (or ``mc``): Monte-Carlo radiation model. This model can be configured to generate macro-photons with :py:data:`radiation_photon_species`.

  This parameter cannot be assigned to photons (mass = 0).

  Radiation is emitted only with the ``"Monte-Carlo"`` model when
  :py:data:`radiation_photon_species` is defined.

.. py:data:: radiation_photon_species

  The :py:data:`name` of the photon species in which the Monte-Carlo :py:data:`radiation_model`
  will generate macro-photons. If unset (or ``None``), no macro-photon will be created.
  The *target* photon species must be have its mass set to 0, and appear *after* the
  particle species in the namelist.

  This parameter cannot be assigned to photons (mass = 0).

.. py:data:: radiation_photon_sampling

  :default: ``1``

  The number of macro-photons generated per emission event, when the macro-photon creation
  is activated (see :py:data:`radiation_photon_species`). The total macro-photon weight
  is still conserved.

  A large number may rapidly slow down the performances and lead to memory saturation.

  This parameter cannot be assigned to photons (mass = 0).
  
.. py:data:: radiation_max_emissions

  :default: ``10``

  The maximum number of emission Monte-Carlo event a macro-particle can undergo during a timestep.
  Since this value is used to allocate some buffers, a high value can saturate memory.

  This parameter cannot be assigned to photons (mass = 0).

.. py:data:: radiation_photon_gamma_threshold

  :default: ``2``

  The threshold on the photon energy for the macro-photon emission when using the
  radiation reaction Monte-Carlo process.
  Under this threshold, the macro-photon from the radiation reaction Monte-Carlo
  process is not created but still taken into account in the energy balance.
  The default value corresponds to twice the electron rest mass energy that
  is the required energy to decay into electron-positron pairs.

  This parameter cannot be assigned to photons (mass = 0).

.. py:data:: relativistic_field_initialization

  :default: ``False``

  Flag for relativistic particles. If ``True``, the electromagnetic fields of this species will added to the electromagnetic fields already present in the simulation.
  This operation will be performed when time equals :py:data:`time_frozen`. See :doc:`/Understand/relativistic_fields_initialization` for details on the computation of the electromagentic fields of a relativistic species.
  To have physically meaningful results, we recommend to place a species which requires this method of field initialization far from other species, otherwise the latter could experience instantly turned-on unphysical forces by the relativistic species' fields.



.. py:data:: multiphoton_Breit_Wheeler

  :default: ``[None,None]``

  An list of the :py:data:`name` of two species: electrons and positrons created through
  the :doc:`/Understand/multiphoton_Breit_Wheeler`.
  By default, the process is not activated.

  This parameter can **only** be assigned to photons species (mass = 0).

.. py:data:: multiphoton_Breit_Wheeler_sampling

  :default: ``[1,1]``

  A list of two integers: the number of electrons and positrons generated per photon decay
  in the :doc:`/Understand/multiphoton_Breit_Wheeler`. The total macro-particle weight is still
  conserved.

  Large numbers may rapidly slow down the performances and lead to memory saturation.

  This parameter can **only** be assigned to photons species (mass = 0).

.. py:data:: keep_interpolated_fields
  
  :default: ``[]``
  
  A list of interpolated fields that should be stored in memory for all particles of this species,
  instead of being located in temporary buffers. These fields can then
  be accessed in some diagnostics such as :ref:`particle binning <DiagParticleBinning>` or
  :ref:`tracking <DiagTrackParticles>`. The available fields are ``"Ex"``, ``"Ey"``, ``"Ez"``, 
  ``"Bx"``, ``"By"`` and ``"Bz"``.
  
  Note that magnetic field components, as they originate from the interpolator,
  are shifted by half a timestep compared to those from the *Fields* diagnostics.
  
  Additionally, the work done by each component of the electric field is available as
  ``"Wx"``, ``"Wy"`` and ``"Wz"``. Contrary to the other interpolated fields, these quantities
  are accumulated over time.

----

.. _Particle_injector:

Particle Injector
^^^^^^^^^^^^^^^^^

Injectors enable to inject macro-particles in the simulation domain from the boundaries.
By default, some parameters that are not specified are inherited from the associated :py:data:`species`.

Each particle injector has to be defined in a ``ParticleInjector`` block::

    ParticleInjector(
        name      = "injector1",
        species   = "electrons1",
        box_side  = "xmin",
        time_envelope = tgaussian(start=0, duration=10., order=4),

        # Parameters inherited from the associated ``species`` by default

        position_initialization = "species",
        momentum_initialization = "rectangular",
        mean_velocity = [0.5,0.,0.],
        temperature = [1e-30],
        number_density = 1,
        particles_per_cell = 16,
    )

.. py:data:: name

    The name you want to give to this injector.
    If you do not specify a name, it will be attributed automatically.
    The name is useful if you want to inject particles at the same position of another injector.

.. py:data:: species

    The name of the species in which to inject the new particles

.. py:data:: box_side

    From where the macro-particles are injected. Options are:

    * ``"xmin"``
    * ``"xmax"``
    * ``"ymin"``
    * ``"ymax"``
    * ``"zmax"``
    * ``"zmin"``

.. py:data:: time_envelope

    :type: a *python* function or a :doc:`time profile <profiles>`
    :default:  ``tconstant()``

    The temporal envelope of the injector.

.. py:data:: position_initialization

    :default: parameters provided the species

    The method for initialization of particle positions. Options are:

    * ``"species"`` or empty ``""``: injector uses the option of the specified :py:data:`species`.
    * ``"regular"`` for regularly spaced. See :py:data:`regular_number`.
    * ``"random"`` for randomly distributed
    * ``"centered"`` for centered in each cell
    * The :py:data:`name` of another injector from which the positions are copied.
      This option requires (1) that the *target* injector's positions are initialized
      using one of the three other options above.

.. py:data:: momentum_initialization

    :default: parameters provided the species

    The method for initialization of particle momenta. Options are:

    * ``"species"`` or empty ``""``: injector uses the option of the specified :py:data:`species`.
    * ``"maxwell-juettner"`` for a relativistic maxwellian (see :doc:`how it is done<maxwell-juttner>`)
    * ``"rectangular"`` for a rectangular distribution

.. py:data:: mean_velocity

    :type: a list of 3 floats or :doc:`profiles <profiles>`
    :default: parameters provided the species

    The initial drift velocity of the particles, in units of the speed of light :math:`c`.

    **WARNING**: For massless particles, this is actually the momentum in units of :math:`m_e c`.

.. py:data:: temperature

    :type: a list of 3 floats or :doc:`profiles <profiles>`
    :default: parameters provided the species

    The initial temperature of the particles, in units of :math:`m_ec^2`.

.. py:data:: particles_per_cell

    :type: float or :doc:`profile <profiles>`
    :default: parameters provided the species

    The number of particles per cell to use for the injector.

.. py:data:: number_density
             charge_density

    :type: float or :doc:`profile <profiles>`
    :default: parameters provided the species

    The absolute value of the number density or charge density (choose one only)
    of the particle distribution, in units of the reference density :math:`N_r` (see :doc:`/Understand/units`)

.. py:data:: regular_number

    :type: A list of as many integers as the simulation dimension

    Same as for :ref:`species`. When ``position_initialization = "regular"``, this sets the number of evenly-spaced
    particles per cell in each direction: ``[Nx, Ny, Nz]`` in cartesian geometries.

----

.. rst-class:: experimental

.. _Particle_merging:

Particle Merging
^^^^^^^^^^^^^^^^

The macro-particle merging method is documented in
the :doc:`corresponding page </Understand/particle_merging>`.
Note that for merging to be able to operate either vectorization or cell sorting must be activated.
It is optionnally specified in the ``Species`` block::

  Species(
      ....

      # Merging
      merging_method = "vranic_spherical",
      merge_every = 5,
      merge_min_particles_per_cell = 16,
      merge_max_packet_size = 4,
      merge_min_packet_size = 4,
      merge_momentum_cell_size = [16,16,16],
      merge_discretization_scale = "linear",
      # Extra parameters for experts:
      merge_min_momentum_cell_length = [1e-10, 1e-10, 1e-10],
      merge_accumulation_correction = True,
  )

.. py:data:: merging_method

  :default: ``"none"``

  The particle merging method to use:

  * ``"none"``: no merging
  * ``"vranic_cartesian"``: method of M. Vranic with a cartesian momentum-space decomposition
  * ``"vranic_spherical"``: method of M. Vranic with a spherical momentum-space decomposition

.. py:data:: merge_every

  :default: ``0``

  Number of timesteps between each merging event
  **or** a :ref:`time selection <TimeSelections>`.

.. py:data:: min_particles_per_cell

  :default: ``4``

  The minimum number of particles per cell for the merging.

.. py:data:: merge_min_packet_size

  :default: ``4``

  The minimum number of particles per packet to merge. Must be greater or equal to 4.

.. py:data:: merge_max_packet_size

  :default: ``4``

  The maximum number of particles per packet to merge.

.. py:data:: merge_momentum_cell_size

  :default: ``[16,16,16]``

  A list of 3 integers defining the number of sub-groups in each direction
  for the momentum-space discretization.

.. py:data:: merge_discretization_scale

  :default: ``"linear"``

  The momentum discretization scale:: ``"linear"`` or ``"log"``.
  The ``"log"`` scale only works with the spherical discretization at the moment.

.. py:data:: merge_min_momentum

  :default: ``1e-5``

  :red:`[for experts]` The minimum momentum value when the log scale
  is chosen (``merge_discretization_scale = log``).
  This avoids a potential 0 value in the log domain.

.. py:data:: merge_min_momentum_cell_length

  :default: ``[1e-10,1e-10,1e-10]``

  :red:`[for experts]` The minimum sub-group length for the momentum-space
  discretization (below which the number of sub-groups is set to 1).

.. py:data:: merge_accumulation_correction

  :default: ``True``

  :red:`[for experts]` Activates the accumulation correction
  (see :doc:`/Understand/particle_merging` for more information).
  The correction only works in linear scale.



----

.. _Lasers:

Lasers
^^^^^^

A laser consists in applying oscillating boundary conditions for the magnetic
field on one of the box sides. The only boundary conditions that support lasers
are ``"silver-muller"`` and ``"PML"`` (see :py:data:`EM_boundary_conditions`).
There are several syntaxes to introduce a laser in :program:`Smilei`:

.. note::

  The following definitions are given for lasers incoming from the ``xmin`` or ``xmax``
  boundaries. For lasers incoming from ``ymin`` or ``ymax``, replace the ``By``
  profiles by ``Bx`` profiles. For lasers incoming from ``zmin`` or ``zmax``,
  replace ``By`` and ``Bz`` profiles by ``Bx`` and ``By`` profiles, respectively.

.. rubric:: 1. Defining a generic wave

..

  .. code-block:: python

    Laser(
        box_side = "xmin",
        space_time_profile = [ By_profile, Bz_profile ]
        space_time_profile_AM = [ Br_mode0, Bt_mode0, Br_mode1, Bt_mode1, ... ]
    )

.. py:data:: box_side

    :default: ``"xmin"``

    Side of the box from which the laser originates: ``"xmin"``, ``"xmax"``, ``"ymin"``,
    ``"ymax"``, ``"zmin"`` or ``"zmax"``.

    In the cases of ``"ymin"`` or ``"ymax"``, replace, in the following profiles,
    coordinates *y* by *x*, and fields :math:`B_y` by :math:`B_x`.

    In the cases of ``"zmin"`` or ``"zmax"``, replace, in the following profiles,
    coordinates *y* by *x*, coordinates *z* by *y*, fields :math:`B_y` by :math:`B_x`
    and fields :math:`B_z` by :math:`B_y`.


.. py:data:: space_time_profile

    :type: A list of two *python* functions

    The full wave expression at the chosen box side. It is a list of **two** *python*
    functions taking several arguments depending on the simulation dimension:
    :math:`(t)` for a 1-D simulation, :math:`(y,t)` for a 2-D simulation (etc.)
    The two functions represent :math:`B_y` and :math:`B_z`, respectively.
    This can be used only in Cartesian geometries.

.. py:data:: space_time_profile_AM

    :type: A list of maximum 2 x ``number_of_AM`` complex valued *python* functions.

    These profiles define the first modes of :math:`B_r` and :math:`B_\theta` of the laser in the
    order shown in the above example. Higher undefined modes are considered zero.
    This can be used only in ``AMcylindrical`` geometry. In this
    geometry a two-dimensional :math:`(x,r)` grid is used and the laser is injected from a
    :math:`x` boundary, thus the provided profiles must be a function of :math:`(r,t)`.



.. rubric:: 2. Defining the wave envelopes

..

  .. code-block:: python

    Laser(
        box_side       = "xmin",
        omega          = 1.,
        chirp_profile  = tconstant(),
        time_envelope  = tgaussian(),
        space_envelope = [ By_profile  , Bz_profile   ],
        phase          = [ PhiY_profile, PhiZ_profile ],
        delay_phase    = [ 0., 0. ]
    )

  This implements a wave of the form:

  .. math::

    B_y(\mathbf{x}, t) = S_y(\mathbf{x})\; T\left(t-t_{0y}\right)
    \;\sin\left( \omega(t) t - \phi_y(\mathbf{x}) \right)

    B_z(\mathbf{x}, t) = S_z(\mathbf{x})\; T\left(t-t_{0z}\right)
    \;\sin\left( \omega(t) t - \phi_z(\mathbf{x}) \right)

  where :math:`T` is the temporal envelope, :math:`S_y` and :math:`S_z` are the
  spatial envelopes, :math:`\omega` is the time-varying frequency,
  :math:`\phi_y` and :math:`\phi_z` are the phases, and we defined the delays
  :math:`t_{0y} = (\phi_y(\mathbf{x})-\varphi_y)/\omega(t)` and
  :math:`t_{0z} = (\phi_z(\mathbf{x})-\varphi_z)/\omega(t)`.

  .. py:data:: omega

    :default: 1.

    The laser angular frequency.

  .. py:data:: chirp_profile

    :type: a *python* function or a :doc:`time profile <profiles>`
    :default: ``tconstant()``

    The variation of the laser frequency over time, such that
    :math:`\omega(t)=` ``omega`` x ``chirp_profile`` :math:`(t)`.

  .. warning::

    This definition of the chirp profile is not standard.
    Indeed, :math:`\omega(t)` as defined here **is not** the instantaneous frequency, :math:`\omega_{\rm inst}(t)`,
    which is obtained from the time derivative of the phase :math:`\omega(t) t`.

    Should one define the chirp as :math:`C(t) = \omega_{\rm inst}(t)/\omega` (with :math:`\omega` defined by the input
    parameter :math:`\mathtt{omega}`), the user can easily obtain the corresponding chirp profile as defined in
    :program:`Smilei` as:

    .. math::

        \mathtt{chirp\_profile}(t) = \frac{1}{t} \int_0^t dt' C(t')\,.

    Let us give as an example the case of a *linear chirp*, with the instantaneous frequency
    :math:`\omega_{\rm inst}(t) = \omega [1+\alpha\,\omega(t-t_0)]`.
    :math:`C(t) = 1+\alpha\,\omega(t-t_0)`. The corresponding input chirp profile reads:

    .. math::

        \mathtt{chirp\_profile}(t) = 1 - \alpha\, \omega t_0 + \frac{\alpha}{2} \omega t

    Similarly, for a *geometric (exponential) chirp* such that :math:`\omega_{\rm inst}(t) = \omega\, \alpha^{\omega t}`,
    :math:`C(t) = \alpha^{\omega t}`, and the corresponding input chirp profile reads:

    .. math::

        \mathtt{chirp\_profile}(t) = \frac{\alpha^{\omega t} - 1}{\omega t \, \ln \alpha}\,.


  .. py:data:: time_envelope

    :type: a *python* function or a :doc:`time profile <profiles>`
    :default:  ``tconstant()``

    The temporal envelope of the laser (field, not intensity).

  .. py:data:: space_envelope

    :type: a list of two *python* functions or two :doc:`spatial profiles <profiles>`
    :default: ``[ 1., 0. ]``

    The two spatial envelopes :math:`S_y` and :math:`S_z`.

  .. py:data:: phase

    :type: a list of two *python* functions or two :doc:`spatial profiles <profiles>`
    :default: ``[ 0., 0. ]``

    The two spatially-varying phases :math:`\phi_y` and :math:`\phi_z`.

  .. py:data:: delay_phase

    :type: a list of two floats
    :default: ``[ 0., 0. ]``

    An extra delay for the time envelopes of :math:`B_y` and :math:`B_z`,
    expressed in terms of phase (:math:`=\omega t`). This delay is applied to the
    :py:data:`time_envelope`, but not to the carrier wave.
    This option is useful in the
    case of elliptical polarization where the two temporal profiles should have a slight
    delay due to the mismatched :py:data:`phase`.



.. rubric:: 3. Defining a 1D planar wave

..

  For one-dimensional simulations, you may use the simplified laser creator::

    LaserPlanar1D(
        box_side         = "xmin",
        a0               = 1.,
        omega            = 1.,
        polarization_phi = 0.,
        ellipticity      = 0.,
        time_envelope    = tconstant(),
        phase_offset     = 0.,
    )

  .. py:data:: a0

    :default: 1.

    The normalized vector potential

  .. py:data:: polarization_phi

    :default: 0.

    The angle of the polarization ellipse major axis relative to the X-Y plane, in radians.

  .. py:data:: ellipticity

    :default: 0.

    The polarization ellipticity: 0 for linear and :math:`\pm 1` for circular.

  .. py:data:: phase_offset
    
    :default: 0.
    
    An extra phase added to both the envelope and to the carrier wave.


.. rubric:: 4. Defining a 2D gaussian wave

..

  For two-dimensional simulations, you may use the simplified laser creator::

    LaserGaussian2D(
        box_side         = "xmin",
        a0               = 1.,
        omega            = 1.,
        focus            = [50., 40.],
        waist            = 3.,
        incidence_angle  = 0.,
        polarization_phi = 0.,
        ellipticity      = 0.,
        time_envelope    = tconstant(),
        phase_offset     = 0.,
    )

  This is similar to ``LaserPlanar1D``, with some additional arguments for
  specific 2D aspects.
  
  .. py:data:: focus

    :type: A list of two floats ``[X, Y]``

    The ``X`` and ``Y`` positions of the laser focus.

  .. py:data:: waist

    The waist value. Transverse coordinate at which the field is at 1/e of its maximum value.

  .. py:data:: incidence_angle

    :default: 0.

    The angle of the laser beam relative to the normal to the injection plane, in radians.


.. rubric:: 5. Defining a 3D gaussian wave

..

  For three-dimensional simulations, you may use the simplified laser creator::

    LaserGaussian3D(
        box_side         = "xmin",
        a0               = 1.,
        omega            = 1.,
        focus            = [50., 40., 40.],
        waist            = 3.,
        incidence_angle  = [0., 0.1],
        polarization_phi = 0.,
        ellipticity      = 0.,
        time_envelope    = tconstant(),
        phase_offset     = 0.,
    )

  This is almost the same as ``LaserGaussian2D``, with the ``focus`` parameter having
  now 3 elements (focus position in 3D), and the ``incidence_angle`` being a list of
  two angles, corresponding to rotations around ``y`` and ``z``, respectively.

  When injecting on ``"ymin"`` or ``"ymax"``, the incidence angles corresponds to
  rotations around ``x`` and ``z``, respectively.

.. rubric:: 6. Defining a gaussian wave with Azimuthal Fourier decomposition

..

  For simulations with ``"AMcylindrical"`` geometry, you may use the simplified laser creator::

    LaserGaussianAM(
        box_side         = "xmin",
        a0               = 1.,
        omega            = 1.,
        focus            = [50.],
        waist            = 3.,
        polarization_phi = 0.,
        ellipticity      = 0.,
        time_envelope    = tconstant()
    )

  Note that here the focus is given in [x] coordinates, since it propagates on the `r=0` axis .

.. rubric:: 7. Defining a generic wave at some distance from the boundary

..

  In some cases, the laser field is not known at the box boundary, but rather at some
  plane inside the box. Smilei can pre-calculate the corresponding wave at the boundary
  using the *angular spectrum method*. This technique is only available in 2D and 3D
  cartesian geometries and requires the python packages *numpy*.
  A :doc:`detailed explanation <laser_offset>` of the method is available.
  The laser is introduced using::

    LaserOffset(
        box_side               = "xmin",
        space_time_profile     = [ By_profile, Bz_profile ],
        offset                 = 10.,
        extra_envelope          = tconstant(),
        keep_n_strongest_modes = 100,
        angle = 10./180.*3.14159
    )

  .. py:data:: space_time_profile

    :type: A list of two *python* functions

    The magnetic field profiles at some arbitrary plane, as a function of space and time.
    The arguments of these profiles are ``(y,t)`` in 2D and ``(y,z,t)`` in 3D.

  .. py:data:: offset

     The distance from the box boundary to the plane where :py:data:`space_time_profile`
     is defined.

  .. py:data:: extra_envelope

    :type: a *python* function or a :doc:`python profile <profiles>`
    :default:  ``lambda *z: 1.``, which means a profile of value 1 everywhere

    An extra envelope applied at the boundary, on top of the :py:data:`space_time_profile`.
    This envelope takes two arguments (``y``, ``t``) in 2D, and three arguments (``y``, ``z``, ``t``)
    in 3D.
    As the wave propagation technique stores a limited number of Fourier modes (in the time
    domain) of the wave, some periodicity can be obtained in the actual laser.
    One may thus observe that the laser pulse is repeated several times.
    The envelope can be used to remove these spurious repetitions.

  .. py:data:: keep_n_strongest_modes

    :default: 100

    The number of temporal Fourier modes that are kept during the pre-processing.
    See :doc:`this page <laser_offset>` for more details.

  .. py:data:: angle

    :default: 0.

    Angle between the boundary and the profile's plane, the rotation being around :math:`z`.
    See :doc:`this page <laser_offset>` for more details.

  .. py:data:: fft_time_window

    :default: :py:data:`simulation_time`

    Time during which the ``space_time_profile`` is sampled (calculating the
    ``LaserOffset`` on the whole simulation duration can be costly). Note that
    the Fourier approach will naturally repeat the signal periodically.
    
  .. py:data:: fft_time_step

    :default: :py:data:`timestep`
    
    Temporal step between each sample of the ``space_time_profile``.
    Chosing a larger step can help reduce the memory load but will remove high temporal frequencies.

  .. py:data:: number_of_processes

    :default: *all available processes*

    The number of MPI processes that will be used for computing the ``LaserOffset``.
    Using more processes computes the FFT faster, but too many processes may
    be very costly in communication. In addition, using too few may not allow
    the arrays to fit in memory.

  .. py:data:: file
  
    :default: ``None``
    
    The path to a ``LaserOffset*.h5`` file generated from a previous simulation. This option
    can help reduce the computation time by re-using the ``LaserOffset`` computation
    from a previous simulation.


----

.. _laser_envelope:

Laser envelope model
^^^^^^^^^^^^^^^^^^^^^^

In all the available geometries, it is possible to model a laser pulse
propagating in the ``x`` direction
using an envelope model (see :doc:`/Understand/laser_envelope` for the advantages
and limits of this approximation).
The fast oscillations of the laser are neglected and all the physical
quantities of the simulation, including the electromagnetic fields and
their source terms, as well as the particles positions and momenta, are
meant as an average over one or more optical cycles.
Effects involving characteristic lengths comparable to the laser central
wavelength (i.e. sharp plasma density profiles) cannot be modeled with
this option.

.. note::

  The envelope model in ``"AMcylindrical"`` geometry is implemented only in the hypothesis of
  cylindrical symmetry, i.e. only one azimuthal mode. Therefore, to use it the user must choose
  ``number_of_AM = 1``.

Contrarily to a standard Laser initialized with the Silver-Müller
boundary conditions, the laser envelope will be entirely initialized inside
the simulation box at the start of the simulation.

Currently only one laser pulse of a given frequency propagating in the positive
`x` direction can be speficified. However, a multi-pulse set-up can be initialized
if a multi-pulse profile is specified, e.g. if the temporal profile is given by two adjacents gaussian functions.
The whole multi-pulse profile would have the same carrier frequency and would propagate in the positive
`x` direction. For the moment it is not possible to specify more than one laser envelope profile, e.g.
two counterpropagating lasers, or two lasers with different carrier frequency.


Please note that describing a laser through its complex envelope loses physical accuracy if its
characteristic space-time variation scales are too small, i.e. of the order of the laser
central wavelength (see :doc:`/Understand/laser_envelope`).
Thus, space-time profiles with variation scales larger than this length should be used.

.. rubric:: 1. Defining a generic laser envelope

..

Following is the generic laser envelope creator ::

    LaserEnvelope(
        omega          = 1.,
        envelope_solver = 'explicit',
        envelope_profile = envelope_profile,
        Envelope_boundary_conditions = [["reflective"]]
        polarization_phi = 0.,
        ellipticity      = 0.
    )


.. py:data:: omega

   :default: ``1.``

   The laser angular frequency.

.. py:data:: envelope_profile

   :type: a *python* function or a :doc:`python profile <profiles>`
   :default: None

   The laser space-time profile, so if the geometry is ``3Dcartesian`` a function of 4 arguments (3 for space, 1 for time) is necessary.
   Please note that the envelope will be entirely initialized in the simulation box
   already at the start of the simulation, so the time coordinate will be applied
   to the ``x`` direction instead of time. It is recommended to initialize the
   laser envelope in vacuum, separated from the plasma, to avoid unphysical
   results.
   Envelopes with variation scales near to the laser wavelength do not
   satisfy the assumptions of the envelope model (see :doc:`/Understand/laser_envelope`),
   yielding inaccurate results.

.. py:data:: envelope_solver

  :default: ``explicit``

  The solver scheme for the envelope equation.

  * ``"explicit"``: an explicit scheme based  on central finite differences.
  * ``"explicit_reduced_dispersion"``: the finite difference derivatives along ``x`` in the ``"explicit"`` solver are substituted by
    optimized derivatives to reduce numerical dispersion. For more accurate results over long distances, the use of this solver is recommended.
    Please note that the CFL limit of this solver is lower than the one of the ``"explicit"`` solver. Thus, a smaller integration 
    timestep may be necessary.

.. py:data:: Envelope_boundary_conditions

  :type: list of lists of strings
  :default: ``[["reflective"]]``

  Defines the boundary conditions used for the envelope. Either ``"reflective"`` or ``"PML"``.
  In the case of ``"PML"``, make sure to define ``"number_of_pml_cells"`` in the ``Main`` block.

.. py:data:: polarization_phi

  :default: 0.

  The angle of the polarization ellipse major axis relative to the X-Y plane, in radians. Needed only for ionization.

.. py:data:: ellipticity

  :default: 0.

  The polarization ellipticity: 0 for linear and 1 for circular. For the moment, only these two polarizations are available.

.. rubric:: 2. Defining a 1D laser envelope

..

Following is the simplified laser envelope creator in 1D ::

    LaserEnvelopePlanar1D(
        a0              = 1.,
        time_envelope   = tgaussian(center=150., fwhm=40.),
        envelope_solver = 'explicit',
        Envelope_boundary_conditions = [ ["reflective"] ],
        polarization_phi = 0.,
        ellipticity      = 0.
    )

.. rubric:: 3. Defining a 2D gaussian laser envelope

..

Following is the simplified gaussian laser envelope creator in 2D ::

    LaserEnvelopeGaussian2D(
        a0              = 1.,
        focus           = [150., 40.],
        waist           = 30.,
        time_envelope   = tgaussian(center=150., fwhm=40.),
        envelope_solver = 'explicit',
        Envelope_boundary_conditions = [ ["reflective"] ],
        polarization_phi = 0.,
        ellipticity      = 0.
    )

.. rubric:: 4. Defining a 3D gaussian laser envelope

..

Following is the simplified laser envelope creator in 3D ::

    LaserEnvelopeGaussian3D(
        a0              = 1.,
        focus           = [150., 40., 40.],
        waist           = 30.,
        time_envelope   = tgaussian(center=150., fwhm=40.),
        envelope_solver = 'explicit',
        Envelope_boundary_conditions = [ ["reflective"] ],
        polarization_phi = 0.,
        ellipticity      = 0.
    )

.. rubric:: 5. Defining a cylindrical gaussian laser envelope

..

Following is the simplified laser envelope creator in ``"AMcylindrical"`` geometry (remember that
in this geometry the envelope model can be used only if ``number_of_AM = 1``) ::

    LaserEnvelopeGaussianAM(
        a0              = 1.,
        focus           = [150.],
        waist           = 30.,
        time_envelope   = tgaussian(center=150., fwhm=40.),
        envelope_solver = 'explicit',
        Envelope_boundary_conditions = [ ["reflective"] ],
        polarization_phi = 0.,
        ellipticity      = 0.
    )


The arguments appearing ``LaserEnvelopePlanar1D``, ``LaserEnvelopeGaussian2D``,
``LaserEnvelopeGaussian3D`` and ``LaserEnvelopeGaussianAM`` have the same meaning they would have in a
normal ``LaserPlanar1D``, ``LaserGaussian2D``, ``LaserGaussian3D`` and ``LaserGaussianAM``,
with some differences:

.. py:data:: time_envelope

   Since the envelope will be entirely initialized in the simulation box
   already at the start of the simulation, the time envelope will be applied
   in the ``x`` direction instead of time. It is recommended to initialize the
   laser envelope in vacuum, separated from the plasma, to avoid unphysical
   results.
   Temporal envelopes with variation scales near to the laser wavelength do not
   satisfy the assumptions of the envelope model (see :doc:`/Understand/laser_envelope`),
   yielding inaccurate results.

.. py:data:: waist

   Please note that a waist size comparable to the laser wavelength does not
   satisfy the assumptions of the envelope model.


It is important to remember that the profile defined through the blocks
``LaserEnvelopePlanar1D``, ``LaserEnvelopeGaussian2D``, ``LaserEnvelopeGaussian3D``
correspond to the complex envelope of the laser vector potential component
:math:`\tilde{A}` in the polarization direction.
The calculation of the correspondent complex envelope for the laser electric field
component in that direction is described in :doc:`/Understand/laser_envelope`.

Note that only order 2 interpolation and projection are supported in presence of
the envelope model for the laser.

The parameters ``polarization_phi`` and ``ellipticity`` specify the polarization state of the laser. In envelope model implemented in :program:`Smilei`,
they are only used to compute the rate of ionization and the initial momentum of the electrons newly created by ionization,
where the polarization of the laser plays an important role (see :doc:`/Understand/ionization`).
For all other purposes (e.g. the particles equations of motions, the computation of the ponderomotive force,
the evolution of the laser), the polarization of the laser plays no role in the envelope model.


----

.. _ExternalField:

External fields
^^^^^^^^^^^^^^^

An initial field can be applied over the whole box
at the beginning of the simulation using the ``ExternalField`` block::

  ExternalField(
      field = "Ex",
      profile = constant(0.01, xvacuum=0.1)
  )

.. py:data:: field

  Field names in Cartesian geometries: ``"Ex"``, ``"Ey"``, ``"Ez"``, ``"Bx"``, ``"By"``, ``"Bz"``, ``"Bx_m"``, ``"By_m"``, ``"Bz_m"``.
  Field names in AM geometry: ``"El_mode_m"``, ``"Er_mode_m"``, ``"Et_mode_m"``, ``"Bl_mode_m"``, ``"Br_mode_m"``, ``"Bt_mode_m"``, ``"Bl_m_mode_m"``, ``"Br_m_mode_m"``, ``"Bt_m_mode_m"``, ``"A_mode_1"``, ``"A0_mode_1"`` .

.. py:data:: profile

  :type: float or :doc:`profile <profiles>`

  The initial spatial profile of the applied field.
  Refer to :doc:`/Understand/units` to understand the units of this field.

  Note that when using standard FDTD schemes, ``B`` fields are given at time ``t=0.5 dt`` and ``B_m`` fields at time ``t=0`` like ``E`` fields.
  It is important to initialize ``B_m`` fields at ``t=0`` if there are particles in the simulation domain at the start of the simulation.
  If ``B_m`` is omited, it is assumed that the magnetic field is constant and that ``B_m=B``.

  Note that in AM geometry all field names must be followed by the number ``"i"`` of the mode that is currently passed with the string ``"_mode_i"``. For instance ``"Er_mode_1"``.
  In this geometry, an external envelope field can also be used. It needs to be initialized at times ``"t=0"`` in ``"A_mode_1"`` and ``"t=-dt"`` in ``"A0_mode_1"``.
  The user must use the ``"_mode_1"`` suffix for these two fields because there is no other possible mode for them.


----

.. _PrescribedField:

Prescribed fields
^^^^^^^^^^^^^^^^^

User-defined electromagnetic fields, with spatio-temporal dependence,
can be superimposed to the self-consistent Maxwell fields.
These fields push the particles but **do not participate in the Maxwell solver**:
they are not self-consistent.
They are however useful to describe charged particles' dynamics in a given
electromagnetic field.

This feature is accessible using the ``PrescribedField`` block::

  from numpy import cos, sin
  def myPrescribedProfile(x,t):
  	return cos(x)*sin(t)

  PrescribedField(
      field = "Ex",
      profile = myPrescribedProfile
  )

.. py:data:: field

  Field names in Cartesian geometries: ``"Ex"``, ``"Ey"``, ``"Ez"``, ``"Bx_m"``, ``"By_m"`` or ``"Bz_m"``.
  Field names in AM geometry: ``"El_mode_m"``, ``"Er_mode_m"``, ``"Et_mode_m"``, ``"Bl_m_mode_m"``, ``"Br_m_mode_m"`` or ``"Bt_m_mode_m"``.

.. warning::

  When prescribing a magnetic field, always use the time-centered fields ``"Bx_m"``, ``"By_m"`` or ``"Bz_m"``.
  These fields are those used in the particle pusher, and are defined at integer time-steps.

.. warning::

  When prescribing a field in AM geometry, the mode "m" must be specified explicitly in the name of the field and the profile
  must return a complex value.

.. warning::

  ``PrescribedFields`` are not visible in the ``Field`` diagnostic, 
  but can be visualised through ``Probes`` and with the fields attributes of ``TrackParticles`` 
  (since they sample the total field acting on the macro-particles).

.. py:data:: profile

  :type: float or :doc:`profile <profiles>`

  The spatio-temporal profile of the applied field: a *python* function
  with arguments (*x*, *t*) or (*x*, *y*, *t*), etc.
  Refer to :doc:`/Understand/units` to understand the units of this field.


----

.. _antennas:

Antennas
^^^^^^^^

An antenna is an extra current applied during the whole simulation.
It is applied using an ``Antenna`` block::

  Antenna(
      field = "Jz",
      space_profile = gaussian(0.01),
      time_profile = tcosine(base=0., duration=1., freq=0.1)
  )

.. py:data:: field

  The name of the current: ``"Jx"``, ``"Jy"`` or ``"Jz"``.

.. py:data:: space_profile

  :type: float or :doc:`profile <profiles>`

  The initial spatial profile of the applied antenna.
  Refer to :doc:`/Understand/units` to understand the units of this current.


.. py:data:: time_profile

  :type: float or :doc:`profile <profiles>`

  The temporal profile of the applied antenna. It multiplies ``space_profile``.

.. py:data:: space_time_profile

  :type: float or :doc:`profile <profiles>`
  
  A space & time profile for the antenna (not compatible with ``space_profile``
  or ``time_profile``). It should have ``N+1``arguments, where ``N`` is the dimension
  of the simulation. For instance ``(x,t)`` in 1D, ``(x,y,t)`` in 2D, etc.
  
  The function must accept ``x``, ``y`` and ``z`` either as floats or numpy arrays.
  If it accepts floats, the return value must be a float.
  If it accepts numpy arrays, these arrays will correspond to the coordinates of 1 patch,
  and the return value must be a numpy array of the same size.



----

Walls
^^^^^

A wall can be introduced using a ``PartWall`` block in order to
reflect, stop, thermalize or kill particles which reach it::

  PartWall(
      kind = "reflective",
      x = 20.
  )

.. py:data:: kind

  The kind of wall: ``"reflective"``, ``"stop"``, ``"thermalize"`` or ``"remove"``.

.. py:data:: x
             y
             z

  Position of the wall in the desired direction. Use only one of ``x``, ``y`` or ``z``.



----

.. _Collisions:

Collisions & reactions
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

:doc:`/Understand/collisions` account for short-range Coulomb interactions of particles (shorter than the 
cell size), but also include other effects such as impact ionization and nuclear reactions.
These are gathered under this section because they are treated as *binary processes* (meaning
they happen during the encounter of two macro-particles).

They are specified by one or several ``Collisions`` blocks::

  Collisions(
      species1 = ["electrons1",  "electrons2"],
      species2 = ["ions1"],
      debug_every = 1000,
      coulomb_log = 0.,
      coulomb_log_factor = 1.,
      ionizing = False,
  #      nuclear_reaction = [],
  )

.. note::

  The screening from bound electrons, which is important when
  the atom is neutral or partially ionized, is accounted for only in the
  case of e-i collisions. To activate it, atom species **must have**
  their :py:data:`atomic_number` defined and non-zero.


.. py:data:: species1
             species2

  Lists of species' :py:data:`name`.

  The collisions and reactions will occur between all species under the group ``species1``
  and all species under the group ``species2``. For example, to collide all
  electrons with ions::

    species1 = ["electrons1", "electrons2"], species2 = ["ions"]

  .. warning::

    This does not make ``electrons1`` collide with ``electrons2``.

  The two groups of species have to be *completely different* OR *exactly equal*.
  In other words, if ``species1`` is not equal to ``species2``,
  then they cannot have any common species.
  If the two groups are exactly equal, we call this situation **intra-collisions**.

  .. note::

    If both lists ``species1`` and ``species2`` contain only one species,
    the algorithm is potentially faster than the situation with several
    species in one or the other list. This is especially true if the
    machine accepts SIMD vectorization.


.. py:data:: every

  :default: 1

  Number of timesteps between each computation of the collisions or reactions.
  Use a number higher than 1 only if you know the collision frequency is low
  with respect to the inverse of the timestep.


.. py:data:: debug_every

  :default: 0

  Number of timesteps between each output of information about collisions or reactions.
  If 0, there will be no outputs.

.. py:data:: time_frozen

  :default: 0.

  The time during which no collisions or reactions happen, in units of :math:`T_r`.
  
.. py:data:: coulomb_log

  :default: 0.

  The Coulomb logarithm.

  * If :math:`= 0`, the Coulomb logarithm is automatically computed for each collision.
  * If :math:`> 0`, the Coulomb logarithm is equal to this value.
  * If :math:`< 0`, collisions are not treated (but other reactions may happen).

.. py:data:: coulomb_log_factor

  :default: 1.

  A constant, strictly positive factor that multiplies the Coulomb logarithm, regardless
  of :py:data:`coulomb_log` being automatically computed or set to a constant value.
  This can help, for example, to compensate artificially-reduced ion masses.

.. _CollisionalIonization:

.. py:data:: ionizing

  :default: ``False``

  :ref:`Collisional ionization <CollIonization>` is set when this parameter is not ``False``.
  It can either be set to the name of a pre-existing electron species (where the ionized
  electrons are created), or to ``True`` (the first electron species in :py:data:`species1`
  or :py:data:`species2` is then chosen for ionized electrons).

  One of the species groups must be all electrons (:py:data:`mass` = 1), and the other
  one all ions of the same :py:data:`atomic_number`.


.. rst-class:: experimental

.. py:data:: nuclear_reaction

  :type: a list of strings
  :default: ``None`` (no nuclear reaction)

  A list of the species names for the products of :ref:`Nuclear reactions <CollNuclearReactions>`
  that may occur during collisions. You may omit product species if they are not necessary
  for the simulation.

  All members of :py:data:`species1` must be the same type of atoms, which is automatically
  recognized by their :py:data:`mass` and :py:data:`atomic_number`. The same applies for
  all members of :py:data:`species2`.

  In the current version, only the reaction D(d,n)He³ is available.

.. rst-class:: experimental

.. py:data:: nuclear_reaction_multiplier

  :type: a float
  :default: 0. (automatically adjusted)

  The rate multiplier for nuclear reactions. It is a positive number that artificially
  increases the occurence of reactions so that a good statistics is obtained. The number
  of actual reaction products is adjusted by changing their weights in order to provide
  a physically correct number of reactions. Leave this number to ``0.`` for an automatic
  rate multiplier: the final number of produced macro-particles will be of the same order
  as that of reactants.



--------------------------------------------------------------------------------

.. _RadiationReaction:

Radiation reaction
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

The block ``RadiationReaction()`` enables to tune the radiation loss properties
(see :doc:`/Understand/radiation_loss`).
Many parameters are used for the generation of the cross-section tables
for the Monte-Carlo emission process.
If the tables already exist in the simulation directory, then they will be read
and no new table will be generated by :program:`Smilei`.
Otherwise, :program:`Smilei` can compute and output these
tables.

::

  RadiationReaction(

    # Radiation parameters
    minimum_chi_continuous = 1e-3,
    minimum_chi_discontinuous = 1e-2,
    table_path = "<path to the external table folder>",

    # Parameters for Niel et al.
    Niel_computation_method = "table",

  )

.. py:data:: minimum_chi_continuous

  :default: 1e-3

  Threshold on the particle quantum parameter *particle_chi*. When a particle has a
  quantum parameter below this threshold, radiation reaction is not taken
  into account.

.. py:data:: minimum_chi_discontinuous

  :default: 1e-2

  Threshold on the particle quantum parameter *particle_chi* between the continuous
  and the discontinuous radiation model.

.. py:data:: table_path

  :default: ``""``

  Path to the **directory** that contains external tables for the radiation losses.
  If empty, the default tables are used.
  Default tables are embedded in the code.
  External tables can be generated using the external tool :program:`smilei_tables` (see :doc:`tables`).

.. py:data:: Niel_computation_method

  :default: ``"table"``

  Method to compute the value of the table *h* of Niel *et al* during the emission process.
  The possible values are:

  * ``"table"``: the *h* function is tabulated. The table is computed at initialization or read from an external file.
  * ``"fit5"``: A polynomial fit of order 5 is used. No table is required.
    The maximal relative error to the reference data is of maximum of 0.02.
    The fit is valid for quantum parameters :math:`\chi` between 1e-3 and 10.
  * ``"fit10"``:  A polynomial fit of order 10 is used. No table is required.
    The precision if better than the fit of order 5 with a maximal relative error of 0.0002.
    The fit is valid for quantum parameters :math:`\chi` between 1e-3 and 10.
  * ``"ridgers"``: The fit of Ridgers given in Ridgers *et al.*, ArXiv 1708.04511 (2017)

  The use of tabulated values is best for accuracy but not for performance.
  Table access prevent total vectorization.
  Fits are vectorizable.

--------------------------------------------------------------------------------

.. _MultiphotonBreitWheeler:

Multiphoton Breit-Wheeler
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

The block ``MultiphotonBreitWheeler`` enables to tune parameters of the
multiphoton Breit-Wheeler process and particularly the table generation.
For more information on this physical mechanism, see :doc:`/Understand/multiphoton_Breit_Wheeler`.

There are three tables used for the multiphoton Breit-Wheeler refers to as the
*integration_dT_dchi*, *min_particle_chi_for_xi* and *xi* table.

::

  MultiphotonBreitWheeler(

    # Path to the tables
    table_path = "<path to the external table folder>",

  )

.. py:data:: table_path

  :default: ``""``

  Path to the **directory** that contains external tables for the multiphoton Breit-Wheeler.
  If empty, the default tables are used.
  Default tables are embedded in the code.
  External tables can be generated using the external tool :program:`smilei_tables` (see :doc:`tables`).

--------------------------------------------------------------------------------

.. _DiagScalar:

*Scalar* diagnostics
^^^^^^^^^^^^^^^^^^^^^

:program:`Smilei` can collect various scalar data, such as total particle energy, total field energy, etc.
This is done by including the block ``DiagScalar``::

  DiagScalar(
      every = 10 ,
      vars = ["Utot", "Ukin", "Uelm"],
      precision = 10
  )

.. py:data:: every

  Number of timesteps between each output **or** a :ref:`time selection <TimeSelections>`.

.. py:data:: vars

  :default: ``[]``

  | List of scalars that will be actually output. Note that most scalars are computed anyways.
  | Omit this argument to include all scalars.

.. py:data:: precision

  :default: 10

  Number of digits of the outputs.

.. warning::

  Scalars diagnostics min/max cell are not yet supported in ``"AMcylindrical"`` geometry.

The full list of available scalars is given in the table below.

.. warning::

  As some of these quantities are integrated in space and/or time, their
  units are unusual, and depend on the simulation dimension.
  All details :ref:`here<integrated_quantities>`.

.. rst-class:: fancy

+--------------------------------------------------------------------------------------------+
| **Space-integrated energy densities**                                                      |
+--------------------------------------------------------------------------------------------+
| +--------------+-------------------------------------------------------------------------+ |
| | Utot         | Total                                                                   | |
| +--------------+-------------------------------------------------------------------------+ |
| | Ukin         | Total kinetic (in the particles)                                        | |
| +--------------+-------------------------------------------------------------------------+ |
| | Uelm         | Total electromagnetic (in the fields)                                   | |
| +--------------+-------------------------------------------------------------------------+ |
| | Uexp         | Expected (Initial :math:`-` lost :math:`+` gained)                      | |
| +--------------+-------------------------------------------------------------------------+ |
| | Ubal         | Balance (Utot :math:`-` Uexp)                                           | |
| +--------------+-------------------------------------------------------------------------+ |
| | Ubal_norm    | Normalized balance (Ubal :math:`/` Utot)                                | |
| +--------------+-------------------------------------------------------------------------+ |
| | Uelm_Ex      | Ex field contribution (:math:`\int E_x^2 dV /2`)                        | |
| +--------------+-------------------------------------------------------------------------+ |
| |              |  ... same for fields Ey, Ez, Bx_m, By_m and Bz_m                        | |
| +--------------+-------------------------------------------------------------------------+ |
| | Urad         | Total radiated                                                          | |
| +--------------+-------------------------------------------------------------------------+ |
| | UmBWpairs    | Total energy converted into electron-position pairs                     | |
| +--------------+-------------------------------------------------------------------------+ |
+--------------------------------------------------------------------------------------------+
| **Space- & time-integrated Energies lost/gained at boundaries**                            |
+--------------------------------------------------------------------------------------------+
| +--------------+-------------------------------------------------------------------------+ |
| | Ukin_bnd     | Time-accumulated kinetic energy exchanged at the boundaries             | |
| +--------------+-------------------------------------------------------------------------+ |
| | Uelm_bnd     | Time-accumulated EM energy exchanged at boundaries                      | |
| +--------------+-------------------------------------------------------------------------+ |
| | PoyXminInst  | Poynting contribution through xmin boundary during the timestep         | |
| +--------------+-------------------------------------------------------------------------+ |
| | PoyXmin      | Time-accumulated Poynting contribution through xmin boundary            | |
| +--------------+-------------------------------------------------------------------------+ |
| |              |  ... same for other boundaries                                          | |
| +--------------+-------------------------------------------------------------------------+ |
| | Ukin_new     | Time-accumulated kinetic energy from new particles (injector)           | |
| +--------------+-------------------------------------------------------------------------+ |
| | Ukin_out_mvw | Time-accumulated kinetic energy lost by the moving window               | |
| +--------------+-------------------------------------------------------------------------+ |
| | Ukin_inj_mvw | Time-accumulated kinetic energy gained by the moving window             | |
| +--------------+-------------------------------------------------------------------------+ |
| | Uelm_out_mvw | Time-accumulated EM energy lost by the moving window                    | |
| +--------------+-------------------------------------------------------------------------+ |
| | Uelm_inj_mvw | Time-accumulated EM energy gained by the moving window                  | |
| +--------------+-------------------------------------------------------------------------+ |
+--------------------------------------------------------------------------------------------+
| **Particle information**                                                                   |
+--------------------------------------------------------------------------------------------+
| +--------------+-------------------------------------------------------------------------+ |
| | Zavg_abc     | Average charge of species "abc" (equals ``nan`` if no particle)         | |
| +--------------+-------------------------------------------------------------------------+ |
| | Dens_abc     |  ... its integrated density                                             | |
| +--------------+-------------------------------------------------------------------------+ |
| | Ukin_abc     |  ... its integrated kinetic energy density                              | |
| +--------------+-------------------------------------------------------------------------+ |
| | Urad_abc     |  ... its integrated radiated energy density                             | |
| +--------------+-------------------------------------------------------------------------+ |
| | Ntot_abc     |  ... and number of macro-particles                                      | |
| +--------------+-------------------------------------------------------------------------+ |
| +--------------+-------------------------------------------------------------------------+ |
+--------------------------------------------------------------------------------------------+
| **Fields information**                                                                     |
+--------------------------------------------------------------------------------------------+
| +--------------+-------------------------------------------------------------------------+ |
| | ExMin        | Minimum of :math:`E_x`                                                  | |
| +--------------+-------------------------------------------------------------------------+ |
| | ExMinCell    |  ... and its location (cell index)                                      | |
| +--------------+-------------------------------------------------------------------------+ |
| | ExMax        | Maximum of :math:`E_x`                                                  | |
| +--------------+-------------------------------------------------------------------------+ |
| | ExMaxCell    |  ... and its location (cell index)                                      | |
| +--------------+-------------------------------------------------------------------------+ |
| |              |  ... same for fields Ey Ez Bx_m By_m Bz_m Jx Jy Jz Rho                  | |
| +--------------+-------------------------------------------------------------------------+ |
+--------------------------------------------------------------------------------------------+

Checkout the :doc:`post-processing <post-processing>` documentation as well.

----

.. _DiagFields:

*Fields* diagnostics
^^^^^^^^^^^^^^^^^^^^

:program:`Smilei` can collect various field data (electromagnetic fields, currents and density)
taken at the location of the PIC grid, both as instantaneous values and averaged values.
This is done by including a block ``DiagFields``::

  DiagFields(
      #name = "my field diag",
      every = 10,
      time_average = 2,
      fields = ["Ex", "Ey", "Ez"],
      #subgrid = None
  )

.. py:data:: name

  Optional name of the diagnostic. Used only for post-processing purposes.

.. py:data:: every

  Number of timesteps between each output **or** a :ref:`time selection <TimeSelections>`.

.. py:data:: flush_every

  :default: 1

  Number of timesteps **or** a :ref:`time selection <TimeSelections>`.

  When ``flush_every`` coincides with ``every``, the output
  file is actually written ("flushed" from the buffer). Flushing
  too often can *dramatically* slow down the simulation.


.. py:data:: time_average

  :default: ``1`` *(no averaging)*

  The number of timesteps for time-averaging.


.. py:data:: fields

  :default: ``[]`` *(all fields are written)*

  List of the field names that are saved. By default, they all are.
  The full list of fields that are saved by this diagnostic:

  .. rst-class:: fancy

  +----------------+-------------------------------------------------------+
  | | Bx           | |                                                     |
  | | By           | | Components of the magnetic field                    |
  | | Bz           | |                                                     |
  +----------------+-------------------------------------------------------+
  | | Bx_m         | |                                                     |
  | | By_m         | | Components of the magnetic field (time-centered)    |
  | | Bz_m         | |                                                     |
  +----------------+-------------------------------------------------------+
  | | Ex           | |                                                     |
  | | Ey           | | Components of the electric field                    |
  | | Ez           | |                                                     |
  +----------------+-------------------------------------------------------+
  | | Jx           | |                                                     |
  | | Jy           | | Components of the total current                     |
  | | Jz           | |                                                     |
  +----------------+-------------------------------------------------------+
  | | Jx_abc       | |                                                     |
  | | Jy_abc       | | Components of the current due to species "abc"      |
  | | Jz_abc       | |                                                     |
  +----------------+-------------------------------------------------------+
  | | Rho          | |  Total charge density                               |
  | | Rho_abc      | |  Charge density of species "abc"                    |
  +----------------+-------------------------------------------------------+

  In ``AMcylindrical`` geometry, the ``x``, ``y`` and ``z``
  indices are replaced by ``l`` (longitudinal), ``r`` (radial) and ``t`` (theta). In addition,
  the angular Fourier modes are denoted by the suffix ``_mode_i`` where ``i``
  is the mode number.
  If a field is specified without its associated mode number, all available modes will be included.
  In summary, the list of fields reads as follows.

  .. rst-class:: fancy

  +------------------------------+-----------------------------------------+
  | | Bl_mode_0, Bl_mode_1, etc. | |                                       |
  | | Br_mode_0, Br_mode_1, etc. | | Components of the magnetic field      |
  | | Bt_mode_0, Bt_mode_1, etc. | |                                       |
  +------------------------------+-----------------------------------------+
  | | El_mode_0, El_mode_1, etc. | |                                       |
  | | Er_mode_0, Er_mode_1, etc. | | Components of the electric field      |
  | | Et_mode_0, Et_mode_1, etc. | |                                       |
  +------------------------------+-----------------------------------------+
  |  The same notation works for Jl, Jr, Jt, and Rho                       |
  +------------------------------+-----------------------------------------+

  In the case of an envelope model for the laser (see :doc:`/Understand/laser_envelope`),
  the following fields are also available:

  .. rst-class:: fancy

  +----------------+-------------------------------------------------------+
  | |              | | Module of laser vector potential's complex envelope |
  | | Env_A_abs    | | :math:`\tilde{A}` (component along the transverse   |
  | |              | | direction)                                          |
  +----------------+-------------------------------------------------------+
  | | Env_Chi      | | Total  susceptibility :math:`\chi`                  |
  +----------------+-------------------------------------------------------+
  | |              | | Module of laser electric field's complex envelope   |
  | | Env_E_abs    | | :math:`\tilde{E}` (component along the transverse   |
  | |              | | direction)                                          |
  +----------------+-------------------------------------------------------+
  | |              | | Module of laser electric field's complex envelope   |
  | | Env_Ex_abs   | | :math:`\tilde{E}_x` (component along the propagation|
  | |              | | direction)                                          |
  +----------------+-------------------------------------------------------+

  In the case the B-TIS3 interpolation is activated (see :doc:`/Understand/algorithms`),
  the following fields are also available:

  .. rst-class:: fancy
       
  +--------------------------------------------+-----------------------------------------------+
  | | By_mBTIS3                                | | Components of the magnetic field            |
  | | By_mBTIS3                                | | for the B-TIS3 interpolation                |
  | |                                          | | (time-centered)                             |
  +--------------------------------------------+-----------------------------------------------+
  | | Br_mBTIS3_mode_0, Br_mBTIS3_mode_1, etc. | | Components of the magnetic field            |
  | | Bt_mBTIS3_mode_0, Bt+mBTIS3_mode_1, etc. | | for the B-TIS3 interpolation                |
  | |                                          | | (``AMcylindrical`` geometry, time-centered) |
  +--------------------------------------------+-----------------------------------------------+


.. Note:: In a given `DiagFields`, all fields must be of the same kind: either real or complex. Therefore To write these last three envelope real fields in ``"AMcylindrical"`` geometry,
          a dedicated block ``DiagFields`` must be defined, e.g. with ``fields = ["Env_A_abs", "Env_Chi"]``.

.. py:data:: subgrid

  :default: ``None`` *(the whole grid is used)*

  A list of slices indicating a portion of the simulation grid to be written by this
  diagnostic. This list must have as many elements as the simulation dimension.
  For example, in a 3D simulation, the list has 3 elements. Each element can be:

  * ``None``, to select the whole grid along that dimension
  * an integer, to select only the corresponding cell index along that dimension
  * a *python* `slice object <https://docs.python.org/3/library/functions.html#slice>`_
    to select regularly-spaced cell indices along that dimension.

  This can be easily implemented using the
  `numpy.s_ expression <https://docs.scipy.org/doc/numpy/reference/generated/numpy.s_.html>`_.
  For instance, in a 3D simulation, the following subgrid selects only every other element
  in each dimension::

    from numpy import s_
    DiagFields( #...
    	subgrid = s_[::2, ::2, ::2]
    )

  while this one selects cell indices included in a contiguous parallelepiped::

    	subgrid = s_[100:300, 300:500, 300:600]


.. py:data:: datatype

  :default: ``"double"``
  
  The data type when written to the HDF5 file. Accepts ``"double"`` (8 bytes) or ``"float"`` (4 bytes).


----

.. _DiagProbe:

*Probe* diagnostics
^^^^^^^^^^^^^^^^^^^

The fields from the previous section are taken at the PIC grid locations,
but it is also possible to obtain the fields at arbitrary locations.
These are called *probes*.

A probe interpolates the fields at either one point (0-D),
several points arranged in a line (1-D),
or several points arranged in a 2-D or 3-D grid.

.. note::

  * **Probes follow the moving window.**
    To obtain the fields at fixed points in the plasma instead, create a cold,
    chargeless species, and :ref:`track the particles <DiagTrackParticles>`.
  * **In "AMcylindrical" geometry**, probes are defined with 3D Cartesian coordinates
    and cannot be separated per mode. Use Field diagnostics for cylindrical coordinates and
    information per mode.
  * **Probes rely on the particle interpolator to compute fields** so that the
    magnetic field is shifted by half a timestep compared to that of *Fields* diagnostics.

To add one probe diagnostic, include the block ``DiagProbe``::

  DiagProbe(
      #name = "my_probe",
      every    = 10,
      origin   = [1., 1.],
      corners  = [
          [1.,10.],
          [10.,1.],
      ],
      number   = [100, 100],
      fields   = ["Ex", "Ey", "Ez"]
  )

.. py:data:: name

  Optional name of the diagnostic. Used only for post-processing purposes.

.. py:data:: every

  Number of timesteps between each output **or** a :ref:`time selection <TimeSelections>`.

.. py:data:: flush_every

  :default: 1

  Number of timesteps **or** a :ref:`time selection <TimeSelections>`.

  When ``flush_every`` coincides with ``every``, the output
  file is actually written ("flushed" from the buffer). Flushing
  too often can *dramatically* slow down the simulation.


.. py:data:: origin

  :type: A list of floats, of length equal to the simulation dimensionality.

  The coordinates of the origin of the probe grid

.. py:data:: corners
             vectors

  :type: A list of lists of floats.

  Defines the corners of the probe grid.
  Each corner is a list of coordinates (as many as the simulation dimensions).

  When using ``corners``, the absolute coordinates of each corner must be specified.
  When using ``vectors``, the coordinates relative to :py:data:`origin` must be specified.

.. py:data:: number

  :type: A list of integers, one for each dimension of the probe.

  The number of points in each probe axis. Must not be defined for a 0-D probe.

.. py:data:: fields

  :default: ``[]``, which means ``["Ex", "Ey", "Ez", "Bx", "By", "Bz", "Jx", "Jy", "Jz", "Rho"]``

  A list of fields among:

  * the electric field components ``"Ex"``, ``"Ey"``, ``"Ez"``
  * the magnetic field components ``"Bx"``, ``"By"``, ``"Bz"``
  * the Poynting vector components ``"PoyX"``, ``"PoyY"``, ``"PoyZ"``
  * the current density components ``"Jx"``, ``"Jy"``, ``"Jz"`` and charge density ``"Rho"``
  * the current density ``"Jx_abc"``, ``"Jy_abc"``, ``"Jz_abc"`` and charge density ``"Rho_abc"``
    of a given species named ``"abc"``

  In the case of an envelope model for the laser (see :doc:`/Understand/laser_envelope`),
  the following fields are also available: ``"Env_Chi"``, ``"Env_A_abs"``, ``"Env_E_abs"``, ``"Env_Ex_abs"``.
  They are respectively the susceptibility, the envelope of the laser transverse vector potential,
  the envelope of the laser transverse electric field and the envelope of the laser longitudinal
  electric field.
  
  If the B-TIS3 interpolation scheme is activated (see :doc:`/Understand/algorithms`),
  the following fields are also available: ``"ByBTIS3"``, ``"BzBTIS3"``.

.. py:data:: time_integral

  :default: ``False``

  If ``True``, the output is integrated over time. As this option forces field interpolation
  at every timestep, it is recommended to use few probe points.

.. py:data:: datatype

  :default: ``"double"``
  
  The data type when written to the HDF5 file. Accepts ``"double"`` (8 bytes) or ``"float"`` (4 bytes).


**Examples of probe diagnostics**

* 0-D probe in 1-D simulation
  ::

    DiagProbe(
        every = 1,
        origin = [1.2]
    )

* 1-D probe in 1-D simulation
  ::

    DiagProbe(
        every = 1,
        origin  = [1.2],
        corners = [[5.6]],
        number  = [100]
    )

* 1-D probe in 2-D simulation
  ::

    DiagProbe(
        every = 1,
        origin  = [1.2, 4.],
        corners = [[5.6, 4.]],
        number  = [100]
    )

* 2-D probe in 2-D simulation
  ::

    DiagProbe(
        every = 1,
        origin   = [0., 0.],
        corners  = [ [10.,0.], [0.,10.] ],
        number   = [100, 100]
    )


----

.. _DiagParticleBinning:

*ParticleBinning* diagnostics
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

A *particle binning diagnostic* collects data from the macro-particles and processes them during runtime.
It does not provide information on individual particles: instead, it produces
**averaged quantities** like the particle density, currents, etc.
The raw data and how it is post-processed by happi is described :doc:`here <binning_units>`.

The data is discretized inside a "grid" chosen by the user. This grid may be of any dimension.

Examples:

* 1-dimensional grid along the position :math:`x` (gives density variation along :math:`x`)
* 2-dimensional grid along positions :math:`x` and :math:`y` (gives density map)
* 1-dimensional grid along the velocity :math:`v_x` (gives the velocity distribution)
* 2-dimensional grid along position :math:`x` and momentum :math:`p_x` (gives the phase-space)
* 1-dimensional grid along the kinetic energy :math:`E_\mathrm{kin}` (gives the energy distribution)
* 3-dimensional grid along :math:`x`, :math:`y` and :math:`E_\mathrm{kin}` (gives the density map for several energies)
* 1-dimensional grid along the charge :math:`Z^\star` (gives the charge distribution)
* 0-dimensional grid (simply gives the total integrated particle density)

Each dimension of the grid is called "axis".

You can add a particle binning diagnostic by including a block ``DiagParticleBinning()`` in the namelist,
for instance::

  DiagParticleBinning(
      #name = "my binning",
      deposited_quantity = "weight",
      every = 5,
      time_average = 1,
      species = ["electrons1", "electrons2"],
      axes = [
          ["x", 0., 10, 100],
          ["ekin", 0.1, 100, 1000, "logscale", "edge_inclusive"]
      ]
  )

.. py:data:: name

  Optional name of the diagnostic. Used only for post-processing purposes.

.. py:data:: deposited_quantity

  The type of data that is summed in each cell of the grid.
  Consider reading :ref:`this <Weights>` to understand the meaning of the ``weight``.

  * ``"weight"`` results in a number density.
  * ``"weight_charge"`` results in a charge density.
  * ``"weight_charge_vx"`` results in the :math:`j_x` current density (same with :math:`y` and :math:`z`).
  * ``"weight_p"`` results in the momentum density (same with :math:`p_x`, :math:`p_y` and :math:`p_z`).
  * ``"weight_ekin"`` results in the energy density.
  * ``"weight_vx_px"`` results in the ``xx`` pressure (same with yy, zz, xy, yz and xz).
  * ``"weight_chi"`` results in the quantum parameter density (only for species with radiation losses).
  * with a user-defined python function, an arbitrary quantity can be calculated (the *numpy*
    module is necessary). This function should take one argument, for instance
    ``particles``, which contains the attributes ``x``, ``y``, ``z``, ``px``, ``py``,
    ``pz``, ``charge``, ``weight``, ``chi`` and ``id`` (additionally, it may also have the
    attributes ``Ex``, ``Bx``, ``Ey``, and so on, depending on :py:data:`keep_interpolated_fields`).
    Each of these attributes is a *numpy* array
    containing the data of all particles in one patch. The function must return a *numpy*
    array of the same shape, containing the desired deposition of each particle. For example,
    defining the following function::

      def stuff(particles):
          return particles.weight * particles.px

    passed as ``deposited_quantity=stuff``, the diagnostic will sum the weights
    :math:`\times\; p_x`.

    You may also pass directly an implicit (*lambda*) function using::

      deposited_quantity = lambda p: p.weight * p.px


.. py:data:: every

  The number of time-steps between each output, **or** a :ref:`time selection <TimeSelections>`.

.. py:data:: flush_every

  :default: 1

  Number of timesteps **or** a :ref:`time selection <TimeSelections>`.

  When ``flush_every`` coincides with ``every``, the output
  file is actually written ("flushed" from the buffer). Flushing
  too often can *dramatically* slow down the simulation.


.. py:data:: time_average

  :default: 1

  The number of time-steps during which the data is averaged. The data is averaged over `time_average` consecutive iterations after the selected time.


.. py:data:: species

  A list of one or several species' :py:data:`name`.
  All these species are combined into the same diagnostic.

.. py:data:: axes

  A list of *axes* that define the grid.
  There may be as many axes as wanted (there may be zero axes).

  Syntax of one axis: ``[type, min, max, nsteps, "logscale", "edge_inclusive"]``

  * ``type`` is one of:

    * ``"x"``, ``"y"``, ``"z"``: spatial coordinates (``"moving_x"`` with a :ref:`moving window<movingWindow>`)
    * ``"px"``, ``"py"``, ``"pz"``, ``"p"``: momenta
    * ``"vx"``, ``"vy"``, ``"vz"``, ``"v"``: velocities
    * ``"gamma"``, ``"ekin"``: energies
    * ``"chi"``: quantum parameter
    * ``"charge"``: the particles' electric charge
    * or a *python function* with the same syntax as the ``deposited_quantity``.
      Namely, this function must accept one argument only, for instance ``particles``,
      which holds the attributes ``x``, ``y``, ``z``, ``px``, ``py``, ``pz``, ``charge``,
      ``weight`` and ``id``. Each of these attributes is a *numpy* array containing the
      data of all particles in one patch. The function must return a *numpy* array of
      the same shape, containing the desired quantity of each particle that will decide
      its location in the histogram binning.

  * The axis is discretized for ``type`` from ``min`` to ``max`` in ``nsteps`` bins.
  * The ``min`` and ``max`` may be set to ``"auto"`` so that they are automatically
    computed from all the particles in the simulation. This option can be bad for performances.
  * The optional keyword ``logscale`` sets the axis scale to logarithmic instead of linear
    (bins become uneven).
  * The optional keyword ``edge_inclusive`` includes the particles outside the range
    [``min``, ``max``] into the extrema bins.

**Examples of particle binning diagnostics**

* Variation of the density of species ``electron1``
  from :math:`x=0` to 1, every 5 time-steps, without time-averaging
  ::

    DiagParticleBinning(
    	deposited_quantity = "weight",
    	every = 5,
    	time_average = 1,
    	species = ["electron1"],
    	axes = [ ["x",    0.,    1.,    30] ]
    )

* Density map from :math:`x=0` to 1, :math:`y=0` to 1
  ::

    DiagParticleBinning(
    	deposited_quantity = "weight",
    	every = 5,
    	time_average = 1,
    	species = ["electron1"],
    	axes = [ ["x",    0.,    1.,    30],
    	         ["y",    0.,    1.,    30] ]
    )

* Velocity distribution from :math:`v_x = -0.1` to :math:`0.1`
  ::

    DiagParticleBinning(
    	deposited_quantity = "weight",
    	every = 5,
    	time_average = 1,
    	species = ["electron1"],
    	axes = [ ["vx",   -0.1,    0.1,    100] ]
    )

* Phase space from :math:`x=0` to 1 and from :math:`px=-1` to 1
  ::

    DiagParticleBinning(
    	deposited_quantity = "weight",
    	every = 5,
    	time_average = 1,
    	species = ["electron1"],
    	axes = [ ["x",    0.,    1.,    30],
    	         ["px",   -1.,   1.,    100] ]
    )

* Energy distribution from 0.01 to 1 MeV in logarithmic scale.
  Note that the input units are :math:`m_ec^2 \sim 0.5` MeV
  ::

    DiagParticleBinning(
    	deposited_quantity = "weight",
    	every = 5,
    	time_average = 1,
    	species = ["electron1"],
    	axes = [ ["ekin",    0.02,    2.,   100, "logscale"] ]
    )

* :math:`x`-:math:`y` density maps for three bands of energy: :math:`[0,1]`, :math:`[1,2]`, :math:`[2,\infty]`.
  Note the use of ``edge_inclusive`` to reach energies up to :math:`\infty`
  ::

    DiagParticleBinning(
    	deposited_quantity = "weight",
    	every = 5,
    	time_average = 1,
    	species = ["electron1"],
    	axes = [ ["x",    0.,    1.,    30],
    	         ["y",    0.,    1.,    30],
    	         ["ekin", 0.,    6.,    3,  "edge_inclusive"] ]
    )

* Charge distribution from :math:`Z^\star =0` to 10
  ::

    DiagParticleBinning(
    	deposited_quantity = "weight",
    	every = 5,
    	time_average = 1,
    	species = ["electron1"],
    	axes = [ ["charge",    -0.5,   10.5,   11] ]
    )


----

.. _DiagScreen:

*Screen* diagnostics
^^^^^^^^^^^^^^^^^^^^

A *screen* collects data from the macro-particles when they cross a surface.
It processes this data similarly to the :ref:`particle binning diagnostics <DiagParticleBinning>`
as it makes a histogram of the macro-particle properties. There are two differences:

* the histogram is made only by the particles that cross the surface
* the data is accumulated for all timesteps.

You can add a screen by including a block ``DiagScreen()`` in the namelist,
for instance::

  DiagScreen(
      #name = "my screen",
      shape = "plane",
      point = [5., 10.],
      vector = [1., 0.],
      direction = "canceling",
      deposited_quantity = "weight",
      species = ["electron"],
      axes = [["a", -10.*l0, 10.*l0, 40],
              ["px", 0., 3., 30]],
      every = 10
  )

.. py:data:: name

  Optional name of the diagnostic. Used only for post-processing purposes.

.. py:data:: shape

   The shape of the screen surface: ``"plane"``, ``"sphere"``, or ``"cylinder"``.

.. py:data:: point

   :type: A list of floats ``[X]`` in 1D,  ``[X,Y]`` in 2D,  ``[X,Y,Z]`` in 3D

   The coordinates of a point that defines the screen surface:
   a point of the ``"plane"``, the center of the ``"sphere"``,
   or a point on the ``"cylinder"`` axis.

.. py:data:: vector

   :type: A list of floats ``[X]`` in 1D,  ``[X,Y]`` in 2D,  ``[X,Y,Z]`` in 3D

   The coordinates of a vector that defines the screen surface:
   the normal to the ``"plane"``, a radius of the ``"sphere"``.
   or the axis of the ``"cylinder"`` (in the latter case, the vector
   norm defines the cylinder radius).

.. py:data:: direction

   :default: ``"both"``

   Determines how particles are counted depending on which side of the screen they come from.

   * ``"both"`` to account for both sides.
   * ``"forward"`` for only the ones in the direction of the ``vector``.
   * ``"backward"`` for only the ones in the opposite direction.
   * ``"canceling"`` to count negatively the ones in the opposite direction.

.. py:data:: deposited_quantity

   Identical to the ``deposited_quantity`` of :ref:`particle binning diagnostics <DiagParticleBinning>`.

.. py:data:: every

  The number of time-steps between each output, **or** a :ref:`time selection <TimeSelections>`.

.. py:data:: flush_every

  :default: 1

  Number of timesteps **or** a :ref:`time selection <TimeSelections>`.

  When ``flush_every`` coincides with ``every``, the output
  file is actually written ("flushed" from the buffer). Flushing
  too often can *dramatically* slow down the simulation.

.. py:data:: species

  A list of one or several species' :py:data:`name`.
  All these species are combined into the same diagnostic.

.. py:data:: axes

  A list of "axes" that define the grid of the histogram.
  It is identical to that of :ref:`particle binning diagnostics <DiagParticleBinning>`, with the
  addition of four types of axes:

  * If ``shape="plane"``, then ``"a"`` and ``"b"`` are the axes perpendicular to the ``vector``.
  * If ``shape="sphere"``, then ``"theta"`` and ``"phi"`` are the angles with respect to the ``vector``.
  * If ``shape="cylinder"``, then ``"a"`` is along the cylinder axis and ``"phi"`` is the angle around it.


----

.. _DiagRadiationSpectrum:

*RadiationSpectrum* diagnostics
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

A *radiation spectrum diagnostic* computes (at a given time) the instantaneous
power spectrum following from the incoherent emission of high-energy
photons by accelerated charge (see :doc:`/Understand/radiation_loss` for more details
on the emission process and its implementation in :program:`Smilei`).

It is similar to the :ref:`particle binning diagnostics <DiagParticleBinning>`,
with an extra axis of binning: the emitted photon energy.
The other axes remain available to the user.

A radiation spectrum diagnostic is defined by a block ``RadiationSpectrum()``::

  DiagRadiationSpectrum(
      #name = "my radiation spectrum",
      every = 5,
      flush_every = 1,
      time_average = 1,
      species = ["electrons1", "electrons2"],
      photon_energy_axis = [0., 1000., 100, 'logscale'],
      axes = []
  )

.. py:data:: name

  Optional name of the diagnostic. Used only for post-processing purposes.

.. py:data:: every

  The number of time-steps between each output, **or** a :ref:`time selection <TimeSelections>`.

.. py:data:: flush_every

  :default: 1

  Number of timesteps **or** a :ref:`time selection <TimeSelections>`.

  When ``flush_every`` coincides with ``every``, the output
  file is actually written ("flushed" from the buffer). Flushing
  too often can *dramatically* slow down the simulation.

.. py:data:: time_average

  :default: 1

  The number of time-steps during which the data is averaged before output.

.. py:data:: species

  A list of one or several species' :py:data:`name` that emit the radiation.
  All these species are combined into the same diagnostic.

.. py:data:: photon_energy_axis

  The axis of photon energies (in units of :math:`m_e c^2`).
  The syntax is similar to that of
  :ref:`particle binning diagnostics <DiagParticleBinning>`.

  Syntax: ``[min, max, nsteps, "logscale"]``

.. py:data:: axes

  An additional list of "axes" that define the grid.
  There may be as many axes as wanted (there may be zero axes).
  Their syntax is the same that for "axes" of a
  :ref:`particle binning diagnostics <DiagParticleBinning>`.


**Examples of radiation spectrum diagnostics**

* Time-integrated over the full duration of the simulation::

    DiagRadiationSpectrum(
        every = Nt,
        time_average = Nt,
        species = ["electrons"],
        photon_energy_axis = [0., 1000., 100, 'logscale'],
        axes = []
    )

* Angularly-resolved instantaneous radiation spectrum.
  The diagnostic considers that all electrons emit radiation in
  the direction of their velocity::

    from numpy import arctan2, pi

    def angle(p):
        return arctan2(p.py,p.px)

    DiagRadiationSpectrum(
        every = 10,
        species = ["electrons"],
        photon_energy_axis = [0., 1000., 100, 'logscale'],
        axes = [
            [angle,-pi,pi,90]
        ]
    )

----

.. _DiagTrackParticles:

*TrackParticles* diagnostics
^^^^^^^^^^^^^^^^^^^^^^^^^^^^

A *particle tracking diagnostic* records the macro-particle positions and momenta at various timesteps.
Typically, this is used for plotting trajectories.

You can add a tracking diagnostic by including a block ``DiagTrackParticles()`` in the namelist,
for instance::

  DiagTrackParticles(
      species = "electron",
      every = 10,
  #    flush_every = 100,
  #    filter = my_filter,
  #    attributes = ["x", "px", "py", "Ex", "Ey", "Bz"]
  )

.. py:data:: species

  The :py:data:`name` of the species to be tracked.

.. py:data:: every

  :default: 0

  Number of timesteps between each output of particles trajectories, **or** a :ref:`time selection <TimeSelections>`.
  If non-zero, the particles positions will be tracked and written in a file named ``TrackParticlesDisordered_abc.h5``
  (where ``abc`` is the species' :py:data:`name`).

.. py:data:: flush_every

  :default: 1

  Number of timesteps **or** a :ref:`time selection <TimeSelections>`.

  When ``flush_every`` coincides with ``every``, the output
  file for tracked particles is actually written ("flushed" from the buffer). Flushing
  too often can *dramatically* slow down the simulation.

.. py:data:: filter

  A python function giving some condition on which particles are tracked.
  If none provided, all particles are tracked.
  To use this option, the `numpy package <http://www.numpy.org/>`_ must
  be available in your python installation.

  The function must have one argument, that you may call, for instance, ``particles``.
  This object has several attributes ``x``, ``y``, ``z``, ``px``, ``py``, ``pz``, ``charge``,
  ``weight`` and ``id`` (additionally, it may also have the
  attributes ``Ex``, ``Bx``, ``Ey``, and so on, depending on :py:data:`keep_interpolated_fields`).
  Each of these attributes
  are provided as **numpy** arrays where each cell corresponds to one particle.
  The function must return a boolean **numpy** array of the same shape, containing ``True``
  for particles that should be tracked, and ``False`` otherwise.

  The following example selects all the particles that verify :math:`-1<p_x<1`
  or :math:`p_z>3`::

    def my_filter(particles):
        return (particles.px>-1.)*(particles.px<1.) + (particles.pz>3.)

.. Note::
  
  * In the ``filter`` function only, the ``px``, ``py`` and ``pz`` quantities
    are not exactly the momenta.
    They are actually the velocities multiplied by the lorentz factor, i.e.,
    :math:`\gamma v_x`, :math:`\gamma v_y` and :math:`\gamma v_z`.
    This is *not* true for the output of the diagnostic.
  * The ``id`` attribute contains the :doc:`particles identification number<ids>`.
    This number is set to 0 at the beginning of the simulation. **Only after particles have
    passed the filter**, they acquire a positive ``id``.
  * For advanced filtration, Smilei provides the quantity ``Main.iteration``,
    accessible within the ``filter`` function. Its value is always equal to the current
    iteration number of the PIC loop. The current time of the simulation is thus
    ``Main.iteration * Main.timestep``.

.. py:data:: attributes

  :default: ``["x","y","z","px","py","pz","w"]``

  A list of strings indicating the particle attributes to be written in the output.
  The attributes may be the particles' spatial coordinates (``"x"``, ``"y"``, ``"z"``),
  their momenta (``"px"``, ``"py"``, ``"pz"``), their electrical charge (``"q"``),
  their statistical weight (``"w"``), their quantum parameter
  (``"chi"``, only for species with radiation losses) or the fields interpolated
  at their  positions (``"Ex"``, ``"Ey"``, ``"Ez"``, ``"Bx"``, ``"By"``, ``"Bz"``).

.. Note:: Here, interpolated fields are normally computed after the Maxwell solver.
  They may thus differ by half a timestep from those computed at the middle of the
  timestep to push particles. When exact values are needed, use the option
  :py:data:`keep_interpolated_fields`.

----

.. rst-class:: experimental

.. _DiagNewParticles:

*NewParticles* diagnostics
^^^^^^^^^^^^^^^^^^^^^^^^^^^^

A *new-particle diagnostic* records the macro-particle information only at the time when
they are generated by :doc:`../Understand/ionization` or other :doc:`../Understand/physics_modules`.

You can add a new-particle diagnostic by including a block ``DiagNewParticles()`` in the namelist,
for instance::

  DiagNewParticles(
      species = "electron",
      every = 10,
  #    attributes = ["x", "px", "py", "Ex", "Ey", "Bz"]
  )

**All the arguments are identical to those of TrackParticles.**
However, there are particular considerations:

* Although the creation of particles is recorded at every timestep, the argument ``every``
  only indicates how often the data is written to the file. It is recommended to avoid 
  small values of ``every`` for better performance.
* In the case of :doc:`../Understand/ionization`, if the chosen ``species`` is that of the ionized electrons,
  then the attribute "``q``" is not the charge of the electron, but the charge of the
  ion, *before ionization occurred*.

----

.. _DiagPerformances:

*Performances* diagnostics
^^^^^^^^^^^^^^^^^^^^^^^^^^^^

The *performances* diagnostic records information on the computational load and timers
for each MPI process  or for each patch in the simulation.

Only one block ``DiagPerformances()`` may be added in the namelist, for instance::

  DiagPerformances(
      every = 100,
  #    flush_every = 100,
  #    patch_information = True,
  )

.. py:data:: every

  :default: 0

  Number of timesteps between each output, **or** a :ref:`time selection <TimeSelections>`.

.. py:data:: flush_every

  :default: 1

  Number of timesteps **or** a :ref:`time selection <TimeSelections>`.

  When ``flush_every`` coincides with ``every``, the output file is actually written
  ("flushed" from the buffer). Flushing too often might *dramatically* slow down the simulation.

.. py:data:: patch_information

  :default: ``False``

  If ``True``, some information is calculated at the patch level (see :py:meth:`Performances`)
  but this may impact the code performances.

----

.. _TimeSelections:

Time selections
^^^^^^^^^^^^^^^

Several components (mainly diagnostics) may require a selection of timesteps to
be chosen by the user. When one of these timesteps is reached, the diagnostics will
output data. A time selection is given through the parameter ``every`` and is a list
of several numbers.

You may chose between five different syntaxes::

  every = [               period                    ] # Syntax 1
  every = [       start,  period                    ] # Syntax 2
  every = [ start,  end,  period                    ] # Syntax 3
  every = [ start,  end,  period,  repeat           ] # Syntax 4
  every = [ start,  end,  period,  repeat,  spacing ] # Syntax 5

where

* ``start`` is the first timestep of the selection (defaults to 0);

* ``end`` is the last timestep of the selection (defaults to ∞);

* ``period`` is the separation between outputs (defaults to 1);

* ``repeat`` indicates how many outputs to do at each period (defaults to 1);

* ``spacing`` is the separation between each repeat (defaults to 1).

For more clarity, this graph illustrates the five syntaxes for time selections:

.. image:: /_static/TimeSelections.png
  :width: 33em
  :align: center

..

.. admonition:: Tips

  * The syntax ``every = period`` is also accepted.
  * Any value set to ``0`` will be replaced by the default value.
  * Special case: ``every=0`` means no output.
  * The numbers may be non-integers (apart from ``repeat``). The closest timesteps are chosen.

----

Profiles
^^^^^^^^^^^^^^^

Some of the quantities described in the previous sections can be profiles that depend on 
space and/or time. See the :doc:`documentation on profiles <profiles>` for detailed
instructions.


----

.. _Checkpoints:

Checkpoints
^^^^^^^^^^^

The simulation state can be saved (*dumped*) at given times (*checkpoints*)
in order to be later *restarted* at that point.

A few things are important to know when you need dumps and restarts.

* Do not restart the simulation in the same directory as the previous one. Files will be
  overwritten, and errors may occur. Create a new directory for your restarted simulation.
* Manage your disk space: each MPI process dumps one file, and the total can be significant.
* The restarted runs must have the same namelist as the initial simulation, except the
  :ref:`Checkpoints` block, which can be modified.

::

  Checkpoints(
      # restart_dir = "dump1",
      dump_step = 10000,
      dump_minutes = 240.,
      exit_after_dump = True,
      keep_n_dumps = 2,
  )

**Parameters to save the state of the current simulation**

  .. py:data:: dump_step

    :default: ``0``

    The number of timesteps between each dump.
    If ``0``, no dump is done.

  .. py:data:: dump_minutes

    :default: ``0.``

    The number of minutes between each dump.
    If ``0.``, no dump is done.

    May be used in combination with :py:data:`dump_step`.

  .. py:data:: exit_after_dump

    :default: ``True``

    If ``True``, the code stops after the first dump. If ``False``, the simulation continues.

  .. py:data:: keep_n_dumps

    :default: ``2``

    This tells :program:`Smilei` to keep, in the current run,  only the last ``n`` dumps.
    Older dumps will be overwritten.

    The default value, ``2``, saves one extra dump in case of a crash during the next dump.

  .. py:data:: file_grouping

    :default: ``0`` (no grouping)

    The maximum number of checkpoint files that can be stored in one directory.
    Subdirectories are created to accomodate for all files.
    This is useful on filesystem with a limited number of files per directory.

  .. py:data:: dump_deflate

    :red:`to do`

**Parameters to restart from a previous simulation**

  .. py:data:: restart_dir

    :default: ``None``

    The directory of a previous run from which :program:`Smilei` should restart.
    For the first run, do not specify this parameter.

    **This path must either absolute or be relative to the current directory.**

    .. Note::

      In many situations, the restarted runs will have the exact same namelist as the initial
      simulation, except this ``restart_dir`` parameter, which points to the previous simulation
      folder.
      You can use the same namelist file, and simply add an extra argument when you launch the
      restart:

      ``mpirun ... ./smilei mynamelist.py "Checkpoints.restart_dir='/path/to/previous/run'"``

  .. py:data:: restart_number

    :default: ``None``

    The number of the dump (in the previous run) that should be used for the restart.
    For the first run, do not specify this parameter.

    In a previous run, the simulation state may have been dumped several times.
    These dumps are numbered 0, 1, 2, etc. until the number :py:data:`keep_n_dumps`.
    In case multiple dumps are kept, the newest one will overwrite the oldest one. 
    To restart the simulation from the most advanced point, specify the dump number 
    corresponding to the newest that was created.


----

Variables defined by Smilei
^^^^^^^^^^^^^^^^^^^^^^^^^^^

:program:`Smilei` passes the following variables to the python interpreter for use in the
namelist. They should not be re-defined by the user!

.. py:data:: smilei_mpi_rank

  The MPI rank of the current process.

.. py:data:: smilei_mpi_size

  The total number of MPI processes.

.. py:data:: smilei_omp_threads

  The number of OpenMP threads per MPI.

.. py:data:: smilei_total_cores

  The total number of cores.

.. note::
  
  These variables can be access during ``happi`` post-processing, e.g.
  ``S.namelist.smilei_mpi_size``.
