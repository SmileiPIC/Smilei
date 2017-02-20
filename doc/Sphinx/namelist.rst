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
        sim_length = [10., 20.], # defines the 2D box dimensions
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

* All quantities are normalized to arbitrary values: see :doc:`units`.

----

Python workflow
^^^^^^^^^^^^^^^

*Python* is started at the beginning of the simulation (one *python* interpreter
for each MPI node). The following steps are executed:

#. A few variables from :program:`Smilei` are passed to *python* so that they are
   available to the user:
   
   * The rank of the current MPI node as :py:data:`smilei_mpi_rank`.
   * The total number of MPI nodes as :py:data:`smilei_mpi_size`.
   * The maximum random integer as :py:data:`smilei_rand_max`.

#. The namelist(s) is executed.

#. *Python* runs :py:data:`cleanup()` if the user has defined it
   (this can be a good place to delete unused heavy variables and unload unused modules).

#. *Python* checks whether the *python* interpreter is needed during the simulation 
   (e.g. the user has defined a temporal :ref:`profile <profiles>` which requires *python*
   to calculate it every timestep). Otherwise, *python* is stopped.

#. If the :py:data:`output_dir` variable was defined, the current working directory
   changes to that value.

All these instructions are summarized in a file ``smilei.py``,
so that the user can directly run ``python -i smilei.py`` for post-processing purposes.

----

Main variables
^^^^^^^^^^^^^^

The block ``Main`` is **mandatory** and has the following syntax::
  
  Main(
      geometry = "1d3v",
      interpolation_order = 2,
      sim_length  = [16. ],
      cell_length = [0.01],
      sim_time    = 15.,
      timestep    = 0.005,
      number_of_patches = [64],
      clrw = 5,
      maxwell_sol = 'Yee',
      bc_em_type_x = ["silver-muller", "silver-muller"],
      bc_em_type_y = ["silver-muller", "silver-muller"],
      time_fields_frozen = 0.,
      referenceAngularFrequency_SI = 0.,
      print_every = 100,
      output_dir = ".",
      random_seed = 0,
  )

.. py:data:: geometry
  
  The geometry of the simulation: ``"1d3v"`` or ``"2d3v"``.
  
  ``1d`` or ``2d`` correspond to the number of spatial dimensions.
  ``3v`` indicates the number of dimensions for velocities.

.. py:data:: interpolation_order
  
  :default: 2
  
  Interpolation order. To this day, only ``2`` is available.


.. py:data:: sim_length
             number_of_cells
  
  A list of floats: size of the simulation box for each dimension of the simulation.
   * Either ``sim_length``, the simulation length in each direction in units of :math:`L_r`,
   * or ``number_of_cells``, the number of cells in each direction.


.. py:data:: cell_length
  
  A list of floats: sizes of one cell in each direction in units of :math:`L_r`.


.. py:data:: sim_time
             number_of_timesteps

  Duration of the simulation.
    * Either ``sim_time``, the simulation duration in units of :math:`T_r`,
    * or ``number_of_timesteps``, the total number of timesteps.


.. py:data:: timestep
             timestep_over_CFL

  Duration of one timestep.
    * Either ``timestep``, in units of :math:`T_r`,
    * or ``timestep_over_CFL``, in units of the *Courant–Friedrichs–Lewy* (CFL) time.


.. py:data:: number_of_patches
  
  A list of integers: the number of patches in each direction.
  Each integer must be a power of 2, and the total number of patches must be
  greater or equal than the number of MPI processes.
  See :doc:`parallelization`.


.. py:data:: clrw
  
  :default: 0.
  
  Cluster width.
  :red:`to do`


.. py:data:: maxwell_sol
  
  :default: 'Yee'
  
  The solver for Maxwell's equations. Only ``"Yee"`` is available at the moment.

.. py:data:: solve_poisson
  
   :default: True
  
   Decides if Poisson correction must be applied or not initially.

.. py:data:: poisson_iter_max
  
  :default: 50000
  
  Maximum number of iteration for the Poisson solver.

.. py:data:: poisson_error_max
  
  :default: 1e-14
  
  Maximum error for the Poisson solver.


.. py:data:: bc_em_type_x
             bc_em_type_y
  
  :type: lists of two strings: ``[bc_min, bc_max]``
  :default: ``["periodic", "periodic"]``
  
  The boundary conditions for the electromagnetic fields.
  The strings ``bc_min`` and ``bc_max`` must be one of the following choices:
  ``"periodic"``, ``"silver-muller"``, or ``"reflective"``.


.. py:data:: time_fields_frozen
  
  :default: 0.
  
  Time, at the beginning of the simulation, during which fields are frozen.


.. _referenceAngularFrequency_SI:

.. py:data:: referenceAngularFrequency_SI
  
  The value of the reference angular frequency :math:`\omega_r` in SI units,
  **only needed when collisions or ionization are requested**.
  This frequency is related to the normalization length according to :math:`L_r\omega_r = c`
  (see :doc:`units`).


.. py:data:: print_every
  
  Number of timesteps between each info output on screen. By default, 10 outputs per
  simulation.


.. py:data:: output_dir

  :default: current working directory
  
  Output directory for the simulation.
  
  **WARNING:** This utility is deprecated and may be removed in a future release.
  Please manage your directories before you run :program:`Smilei`.

.. py:data:: random_seed

  :default: the machine clock

  The value of the random seed. To create a per-processor random seed, you may use
  the variable  :py:data:`smilei_mpi_rank`.

----

Load Balancing
^^^^^^^^^^^^^^

The block ``LoadBalancing`` is optional. If you do not define it, load balancing will
occur every 150 iterations.

.. code-block:: python
  
  LoadBalancing(
      initial_balance = True
      every = 150,
      coef_cell = 1.,
      coef_frozen = 0.1,
  )

.. py:data:: initial_balance
  
  :default: True
  
  Decides if the load must be balanced at initialization. If not, the same amount of
  patches will be attributed to each MPI rank.

.. py:data:: every
  
  :default: 150
  
  An integer: the number of timesteps between each load balancing (patches are
  exchanged between MPI processes to reduce load imbalance).
  
.. py:data:: coef_cell
  
  :default: 1.
  
  :red:`to do`
  
.. py:data:: coef_frozen
  
  :default: 0.1
  
  :red:`to do`


----

Moving window
^^^^^^^^^^^^^

The block ``MovingWindow`` is optional. The window does not move it you do not define it.

.. code-block:: python
  
  MovingWindow(
      time_start = 0.,
      velocity_x = 1.,
  )


.. py:data:: time_start

  :default: 0.
  
  The time at which the window starts moving.


.. py:data:: velocity_x

  :default: 0.
  
  The velocity of the moving window in the `x` direction.

----

.. _Species:

Species
^^^^^^^

Each species has to be defined in a ``Species`` block::

  Species(
      species_type      = "electrons1",
      initPosition_type = "random",
      initMomentum_type = "maxwell-juettner",
      n_part_per_cell = 100,
      mass = 1.,
      atomic_number = None,
      nb_density = 10.,
      # charge_density = None,
      charge = -1.,
      mean_velocity = [0.],
      temperature = [1e-10],
      bc_part_type_xmin = "refl",
      bc_part_type_xmax = "refl",
      # bc_part_type_ymax = None,
      # bc_part_type_ymin = None,
      # thermT = None,
      # thermVelocity = None,
      time_frozen = 0.0,
      # ionization_model = "none",
      # ionization_electrons = None,
      # radiating = False,
      isTest = False,
      track_every = 10,
      track_ordered = False,
      track_flush_every = 100,
      c_part_max = 1.0,
      dynamics_type = "norm",
  )

.. py:data:: species_type
  
  The name you want to give to this species.

.. py:data:: initPosition_type
  
   The initialization of particle positions:
   
   * ``"regular"`` for regularly spaced
   * ``"random"`` for randomly distributed
   * ``"centered"`` for centered in each cell


.. py:data:: initMomentum_type
  
  The initialization of particle momenta:
  
  * ``"maxwell-juettner"`` for a relativistic maxwellian (see :doc:`how it is done<maxwell-juttner>`)
  * ``"rectangular"`` for a rectangular distribution
  * ``"cold"`` for zero temperature
  
  The first 2 distributions depend on the parameter :py:data:`temperature` explained below.

.. py:data:: n_part_per_cell
  
  :type: float or *python* function (see section :ref:`profiles`)
  
  The number of particles per cell.


.. py:data:: mass
  
  The mass of particles, in units of the electron mass :math:`m_e`.


.. py:data:: atomic_number
  
  :default: 0

  The atomic number of the particles, required only for ionization.
  It must be lower than 101.


.. py:data:: nb_density
             charge_density
  
  :type: float or *python* function (see section :ref:`profiles`)
  
  The absolute value of the number density or charge density (choose one only)
  of the particle distribution, in units of the reference density :math:`N_r` (see :doc:`units`).


.. py:data:: charge
  
  :type: float or *python* function (see section :ref:`profiles`)
  
  The particle charge, in units of the electron charge :math:`e`.


.. py:data:: mean_velocity
  
  :type: a list of 3 floats or *python* functions (see section :ref:`profiles`)
  
  The initial drift velocity of the particles, in units of the speed of light :math:`c`.


.. py:data:: temperature
  
  :type: a list of 3 floats or *python* functions (see section :ref:`profiles`)
  
  The initial temperature of the particles, in units of :math:`m_ec^2`.


.. py:data:: bc_part_type_xmin
             bc_part_type_xmax
             bc_part_type_ymin
             bc_part_type_ymax
  
  The boundary condition for particles: ``"refl"`` for *reflecting*, ``"supp"`` for
  *suppressing*, ``"stop"`` for *stopping*, ``"periodic"``, and ``"thermalize"``.
  
.. py:data:: thermT
  
  :default: None
  
  :red:`to do`

.. py:data:: thermVelocity
  
  :default: None
  
  :red:`to do`

.. py:data:: time_frozen
  
  :default: 0.
  
  The time during which the particle positions are not updated, in units of :math:`T_r`.


.. py:data:: ionization_model
  
  :default: ``"none"``
  
  The model for field ionization. Currently, only ``"tunnel"`` is available.
  See :ref:`this <CollisionalIonization>` for collisional ionization instead.


.. py:data:: ionization_electrons
  
  The name of the electron species that field ionization uses when creating new electrons.


.. py:data:: radiating
  
  :default: ``False``
  
  :red:`to do`


.. py:data:: isTest
  
  :default: ``False``
  
  Flag for test particles. If ``True``, this species will contain only test particles
  which do not participate in the charge and currents.

.. py:data:: track_every
  
  :default: 0
  
  Number of timesteps between each output of particles trajectories, **or** a :ref:`time selection <TimeSelections>`.
  If non-zero, the particles positions will be tracked and written in a file named ``TrackParticles_abc.h5``
  (where ``abc`` is :py:data:`species_type`).

.. py:data:: track_ordered
  
  :default: False
  
  If ``True``, tracked particles will be sorted by their ID at run-time. This can be very
  slow. If ``False``, the sorting will occur at post-processing. Again, this may be slow,
  but better for most simulations.

.. py:data:: track_flush_every
  
  :default: 1
  
  Number of timesteps **or** a :ref:`time selection <TimeSelections>`.
  
  When :py:data:`track_flush_every` coincides with :py:data:`track_every`, the output
  file for tracked particles is actually written ("flushed" from the buffer). Flushing
  too often can *dramatically* slow down the simulation.

.. py:data:: c_part_max
  
  :red:`to do`


.. py:data:: dynamics_type
  
  :red:`to do`



----

Lasers
^^^^^^

A laser consists in applying oscillating boundary conditions for the magnetic
field on one of the box sides. The only boundary conditions that support lasers
are ``"silver-muller"`` (see :py:data:`bc_em_type_x`).
There are several syntaxes to introduce a laser in :program:`Smilei`:

.. rubric:: 1. Defining a generic wave

..

  .. code-block:: python
    
    Laser(
        boxSide = "xmin",
        space_time_profile = [ By_profile, Bz_profile ]
    )
  
  .. py:data:: boxSide
    
    :default: ``"xmin"``
    
    Side of the box from which the laser originates: at the moment, only ``"xmin"`` and
    ``"xmax"`` are supported.
    
  .. py:data:: space_time_profile
  
    :type: A list of two *python* functions
    
    The full wave expression at the chosen box side. It is a list of **two** *python*
    functions taking several arguments depending on the simulation dimension:
    :math:`(t)` for a 1-D simulation, :math:`(y,t)` for a 2-D simulation (etc.)
    The two functions represent :math:`B_y` and :math:`B_z`, respectively.
  

.. rubric:: 2. Defining the wave envelopes

..
  
  .. code-block:: python
    
    Laser(
        boxSide        = "xmin",
        omega          = 1.,
        chirp_profile  = tconstant(),
        time_envelope  = tgaussian(),
        space_envelope = [ By_profile  , Bz_profile   ],
        phase          = [ PhiY_profile, PhiZ_profile ]
    )
  
  This implements a wave of the form:
  
  .. math::
    
    B_y(\mathbf{x}, t) = S_y(\mathbf{x})\; T\left[t-\phi_y(\mathbf{x})/\omega(t)\right]
    \;\sin\left( \omega(t) t - \phi_y(\mathbf{x}) \right)
    
    B_z(\mathbf{x}, t) = S_z(\mathbf{x})\; T\left[t-\phi_z(\mathbf{x})/\omega(t)\right]
    \;\sin\left( \omega(t) t - \phi_z(\mathbf{x}) \right)
  
  where :math:`T` is the temporal envelope, :math:`S_y` and :math:`S_y` are the
  spatial envelopes, :math:`\omega` is the time-varying frequency, and 
  :math:`\phi_y` and :math:`\phi_z` are the phases.
  
  .. py:data:: omega
    
    :default: 1.
    
    The laser angular frequency.
    
  .. py:data:: chirp_profile
    
    :type: a *python* function or a :ref:`time profile <profiles>`
    :default: ``tconstant()``
    
    The variation of the laser frequency over time, such that
    :math:`\omega(t)=\mathtt{omega}\times\mathtt{chirp\_profile}(t)`.
    
  .. py:data:: time_envelope
    
    :type: a *python* function or a :ref:`time profile <profiles>`
    :default:  ``tconstant()``
    
    The temporal envelope of the laser.
    
  .. py:data:: space_envelope
    
    :type: a list of two *python* functions or two :ref:`spatial profiles <profiles>`
    :default: ``[ 1., 0. ]``
    
    The two spatial envelopes :math:`S_y` and :math:`S_z`.
    
  .. py:data:: phase
    
    :type: a list of two *python* functions or two :ref:`spatial profiles <profiles>`
    :default: ``[ 0., 0. ]``
    
    The two spatially-varying phases :math:`\phi_y` and :math:`\phi_z`.



.. rubric:: 3. Defining a 1D planar wave

..

  For one-dimensional simulations, you may use the simplified laser creator::
    
    LaserPlanar1D(
        boxSide         = "xmin",
        a0              = 1.,
        omega           = 1.,
        polarizationPhi = 0.,
        ellipticity     = 0.,
        time_envelope   = tconstant()
    )
  
  .. py:data:: a0
  
    :default: 1.
    
    The normalized vector potential
    
  .. py:data:: polarizationPhi
    
    :default: 0.
    
    The angle of the polarization ellipse major axis relative to the X-Y plane, in radians.
    
  .. py:data:: ellipticity
    
    :default: 0.
    
    The polarization ellipticity: 0 for linear and :math:`\pm 1` for circular.



.. rubric:: 4. Defining a 2D gaussian wave

..

  For two-dimensional simulations, you may use the simplified laser creator::
    
    LaserGaussian2D(
        boxSide         = "xmin",
        a0              = 1.,
        omega           = 1.,
        focus           = [50., 40.],
        waist           = 3.,
        incidence_angle = 0.,
        polarizationPhi = 0.,
        ellipticity     = 0.,
        time_envelope   = tconstant()
    )
  
  .. py:data:: focus
    
    :type: A list of two floats ``[X, Y]``
    
    The ``X`` and ``Y`` positions of the laser focus.
    
  .. py:data:: waist
    
    The waist value. Transverse coordinate at which the field is at 1/e of its maximum value.
    
  .. py:data:: incidence_angle
    
    :default: 0.
    
    The angle of the laser beam relative to the X axis, in radians.
  
  .. py:data:: time_envelope
    
     Time envelope of the field (not intensity).


.. rubric:: 5. Defining a 3D gaussian wave

..

  For three-dimensional simulations, you may use the simplified laser creator::
    
    LaserGaussian3D(
        boxSide         = "xmin",
        a0              = 1.,
        omega           = 1.,
        focus           = [50., 40., 40.],
        waist           = 3.,
        incidence_angle = [0., 0.1], 
        polarizationPhi = 0.,
        ellipticity     = 0.,
        time_envelope   = tconstant()
    )
  
  This is almost the same as ``LaserGaussian2D``, with the ``focus`` parameter having
  now 3 elements (focus position in 3D), and the ``incidence_angle`` being a list of
  two angles, corresponding to rotations around `y` and `z`, respectively.



----

.. _ExtField:

External fields
^^^^^^^^^^^^^^^

An external field can be applied using an ``ExtField`` block::

  ExtField(
      fields = ["Ex"],
      profile = constant(0.01, xvacuum=0.1)
  )

.. py:data:: field
               
  Field name: ``"Ex"``, ``"Ey"``, ``"Ez"``, ``"Bx"``, ``"By"`` or ``"Bz"``.
  
.. py:data:: profile
  
  :type: float or *python* function (see section :ref:`profiles`)
  
  The initial spatial profile of the applied field.
  Refer to :doc:`units` to understand the units of this field.


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
  
  :type: float or *python* function (see section :ref:`profiles`)
  
  The initial spatial profile of the applied antenna.
  Refer to :doc:`units` to understand the units of this current.


.. py:data:: time_profile
  
  :type: float or *python* function (see section :ref:`profiles`)
  
  The temporal profile of the applied antenna. It multiplies ``space_profile``.


----

.. _profiles:

Profiles
^^^^^^^^

Several quantities require the input of a profile: particle charge, particle density,
external fields, etc. Depending on the case, they can be *spatial* or *temporal*
profiles.

.. rubric:: 1. Constant profiles

* ``Species( ... , charge = -3., ... )`` defines a species with charge :math:`Z^\star=3`.

* ``Species( ... , nb_density = 10., ... )`` defines a species with density :math:`10\,N_r`.
  You can choose ``nb_density`` (*number density*) or ``charge_density``

* ``Species( ... , mean_velocity = [0.05, 0., 0.], ... )`` defines a species
  with drift velocity :math:`v_x = 0.05\,c` over the whole box.

* ``Species(..., initMomentum_type="maxwell-juettner", temperature=[1e-5], ...)`` defines
  a species with a Maxwell-Jüttner distribution of temperature :math:`T = 10^{-5}\,m_ec^2` over the whole box.
  Note that the temperature may be anisotropic: ``temperature=[1e-5, 2e-5, 2e-5]``.

* ``Species( ... , n_part_per_cell = 10., ... )`` defines a species with 10 particles per cell.

* ``ExtField( field="Bx", profile=0.1 )`` defines a constant external field :math:`B_x = 0.1 B_r`.


.. rubric:: 2. *Python* profiles

..

  Any *python* function can be a profile. You must have basic *python* knowledge to build these functions.
  
  Examples::
  
    def f(x):
        if x<1.: return 0.
        else: return 1.
  
  .. code-block:: python
  
    def f(x,y):    # two variables for 2D simulation
        import math
        twoPI = 2.* math.pi
        return math.cos(  twoPI * x/3.2 )
  
  .. code-block:: python
    
    f = lambda x: x**2 - 1
  
  
  
  Once the function is created, you have to include it in the block you want,
  for example::
  
    Species( ... , charge = f, ... )
    
    Species( ... , mean_velocity = [f, 0, 0], ... )
  

.. rubric:: 3. Pre-defined *spatial* profiles

..

  .. py:function:: constant(value, xvacuum=0., yvacuum=0.)
  
    :param value: the magnitude
    :param xvacuum: vacuum region before the start of the profile.
  
  .. py:function:: trapezoidal(max, \
            xvacuum=0., xplateau=None, xslope1=0., xslope2=0., \
            yvacuum=0., yplateau=None, yslope1=0., yslope2=0. )
  
    :param max: maximum value
    :param xvacuum: empty length before the ramp up
    :param xplateau: length of the plateau (default is :py:data:`sim_length` :math:`-` ``xvacuum``)
    :param xslope1: length of the ramp up
    :param xslope2: length of the ramp down
  
  .. py:function:: gaussian(max, \
     xvacuum=0., xlength=None, xfwhm=None, xcenter=None, xorder=2, \
     yvacuum=0., ylength=None, yfwhm=None, ycenter=None, yorder=2 )
  
    :param max: maximum value
    :param xvacuum: empty length before starting the profile
    :param xlength:  length of the profile (default is :py:data:`sim_length` :math:`-` ``xvacuum``)
    :param xfwhm: gaussian FWHM (default is ``xlength/3.``)
    :param xcenter: gaussian center position (default is in the middle of ``xlength``)
    :param xorder: order of the gaussian.
    :note: If ``yorder`` equals 0, then the profile is constant over :math:`y`.
  
  .. py:function:: polygonal( xpoints=[], xvalues=[] )
  
    :param xpoints: list of the positions of the points
    :param xvalues: list of the values of the profile at each point
  
  .. py:function:: cosine( base, amplitude=1., \
           xvacuum=0., xlength=None, xphi=0., xnumber=1 )
  
    :param base: offset of the profile value
    :param amplitude: amplitude of the cosine
    :param xvacuum: empty length before starting the profile
    :param xlength: length of the profile (default is :py:data:`sim_length` :math:`-` ``xvacuum``)
    :param xphi: phase offset
    :param xnumber: number of periods within ``xlength``
  
  .. py:function:: polynomial( x0=0., y0=0., z0=0., order0=[], order1=[], ... )
    
    :param x0,y0: The reference position(s)
    :param order0: Coefficient for the 0th order
    :param order1: Coefficient for the 1st order (2 coefficients in 2D)
    :param order2: Coefficient for the 2nd order (3 coefficients in 2D)
    :param etc:
    
    Creates a polynomial of the form
    
    .. math::
      
      \begin{eqnarray}
      &\sum_i a_i(x-x_0)^i & \quad\mathrm{in\, 1D}\\
      &\sum_i \sum_j a_{ij}(x-x0)^{i-j}(y-y0)^j & \quad\mathrm{in\, 2D}\\
      &\sum_i \sum_j \sum_k a_{ijk}(x-x0)^{i-j-k}(y-y0)^j(z-z0)^k & \quad\mathrm{in\, 3D}
      \end{eqnarray}
    
    Each ``orderi`` is a coefficient (or list of coefficents) associated to the order ``i``.
    In 1D, there is only one coefficient per order. In 2D, each ``orderi`` is a list
    of ``i+1`` coefficients. For instance, the second order has three coefficients
    associated to :math:`x^2`, :math:`xy` and :math:`y^2`, respectively.
    In 3D, each ``orderi`` is a list of ``(i+1)*(i+2)/2`` coefficients. For instance,
    the second order has 6 coefficients associated to :math:`x^2`, :math:`xy`, :math:`xz`,
    :math:`y^2`, :math:`yz` and :math:`z^2`, respectively.
  
  **Examples**::
    
    Species( ... , density = gaussian(10., xfwhm=0.3, xcenter=0.8), ... )
    
    ExtField( ..., profile = constant(2.2), ... )


.. rubric:: 4. Pre-defined *temporal* profiles

..

  .. py:function:: tconstant(start=0.)
  
    :param start: starting time
  
  .. py:function:: ttrapezoidal(start=0., plateau=None, slope1=0., slope2=0.)
  
    :param start: starting time
    :param plateau: duration of the plateau (default is :py:data:`sim_time` :math:`-` ``start``)
    :param slope1: duration of the ramp up
    :param slope2: duration of the ramp down
  
  .. py:function:: tgaussian(start=0., duration=None, fwhm=None, center=None, order=2)
  
    :param start: starting time
    :param duration: duration of the profile (default is :py:data:`sim_time` :math:`-` ``start``)
    :param fwhm: gaussian FWHM (default is ``duration/3.``)
    :param center: gaussian center time (default is in the middle of ``duration``)
    :param order: order of the gaussian
  
  .. py:function:: tpolygonal( points=[], values=[] )
  
    :param points: list of times
    :param values: list of the values at each time
  
  .. py:function:: tcosine( base=0., amplitude=1., start=0., duration=None, phi=0., freq=1. )
  
    :param base: offset of the profile value
    :param amplitude: amplitude of the cosine
    :param start: starting time
    :param duration: duration of the profile (default is :py:data:`sim_time` :math:`-` ``start``)
    :param phi: phase offset
    :param freq: frequency
  
  .. py:function:: tpolynomial( t0=0., order0=[], order1=[], ... )
    
    :param t0: The reference position
    :param order0: Coefficient for the 0th order
    :param order1: Coefficient for the 1st order
    :param order2: Coefficient for the 2nd order
    :param etc:
    
    Creates a polynomial of the form :math:`\sum_i a_i(t-t_0)^i`.
  
  **Example**::
    
    Antenna( ... , time_profile = tcosine(freq=0.01), ... )


.. rubric:: Illustrations of the pre-defined spatial and temporal profiles
  
.. image:: _static/pythonprofiles.png

.. image:: _static/pythonprofiles_t.png


----

Walls
^^^^^

A wall can be introduced using a ``PartWall`` block in order to
reflect, stop, thermalize or kill particles which reach it::

  PartWall(
      kind = "refl",
      x = 20.
  )

.. py:data:: kind
  
  The kind of wall: ``"refl"``, ``"stop"``, ``"thermalize"`` or ``"supp"``;
  corresponding to a *reflective*, *stopping*, *thermalizing* or *suppressing* wall,
  respectively.
  
.. py:data:: x
             y
             z
  
  Position of the wall in the desired direction. Use only one of ``x``, ``y`` or ``z``.



----

.. _Collisions:

Collisions
^^^^^^^^^^

To have binary collisions in :program:`Smilei`, add one or several ``Collisions`` blocks::

  Collisions(
      species1 = ["electrons1",  "electrons2"],
      species2 = ["ions1"],
      coulomb_log = 5.,
      debug_every = 1000,
      ionizing = False,
  )


.. py:data:: species1
             species2
  
  Lists of species names (see :py:data:`species_type`).
  
  The collisions will occur between all species under the group ``species1``
  and all species under the group ``species2``. For example, to collide all
  electrons with ions::
    
    species1 = ["electrons1", "electrons2"], species2 = ["ions"]

  .. warning::
    
    This does not make ``electrons1`` collide with ``electrons2``.
  
  The two groups of species have to be *completely different* OR *exactly equal*.
  In other words, if ``species1`` is not equal to ``species2``,
  then they cannot have any common species.
  If the two groups are exactly equal, we call this situation **intra-collisions**.


.. py:data:: coulomb_log
  
  :default: 0.
  
  The Coulomb logarithm.
  
  * If :math:`= 0`, the Coulomb logarithm is automatically computed for each collision.
  * If :math:`> 0`, the Coulomb logarithm is equal to this value.


.. py:data:: debug_every
  
  :default: 0
  
  | Number of timesteps between each output of information about collisions.
  | If 0, there will be no outputs.


.. _CollisionalIonization:

.. py:data:: ionizing
  
  :default: False
  
  If ``True``, :ref:`collisional ionization <CollIonization>` will occur. One of the 
  species groups must be all electrons (:py:data:`mass` = 1), and the other one all ions of the
  same :py:data:`atomic_number`.


For more details about the collision scheme in :program:`Smilei`, see :doc:`collisions`


----

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



The full list of scalars that are saved by this diagnostic:


.. rst-class:: nowrap

+----------------+---------------------------------------------------------------------------+
| **Global energies**                                                                        |
+----------------+---------------------------------------------------------------------------+
| | Utot         | | Total energy                                                            |
| | Ukin         | | Total kinetic energy (in the particles)                                 |
| | Uelm         | | Total EM energy (in the fields)                                         |
| | Uexp         | | Expected value (Initial energy :math:`-` lost :math:`+` gained)         |
| | Ubal         | | Energy balance (Utot :math:`-` Uexp)                                    |
| | Ubal_norm    | | Normalized energy balance (Ubal :math:`/` Utot)                         |
| | Uelm_Ex      | | Energy in Ex field (:math:`\int E_x^2 dV /2`)                           |
| |              | |  ... and idem for fields Ey, Ez, Bx_m, By_m and Bz_m                    |
+----------------+---------------------------------------------------------------------------+
| **Energies lost/gained at boundaries**                                                     |
+----------------+---------------------------------------------------------------------------+
| | Ukin_bnd     | | Kinetic energy exchanged at the boundaries during the timestep          |
| | Uelm_bnd     | | EM energy exchanged at boundaries during the timestep                   |
| | Ukin_out_mvw | | Kinetic energy lost during the timestep due to the moving window        |
| | Ukin_inj_mvw | | Kinetic energy injected during the timestep due to the moving window    |
| | Uelm_out_mvw | | EM energy lost during the timestep due to the moving window             |
| | Uelm_inj_mvw | | EM energy injected during the timestep due to the moving window         |
+----------------+---------------------------------------------------------------------------+
| **Species information**                                                                    |
+----------------+---------------------------------------------------------------------------+
| | Dens_abc     | | Average density of species "abc"                                        |
| | Zavg_abc     | |  ... its average charge                                                 |
| | Ukin_abc     | |  ... its total kinetic energy                                           |
| | Ntot_abc     | |  ... and number of particles                                            |
+----------------+---------------------------------------------------------------------------+
| **Fields information**                                                                     |
+----------------+---------------------------------------------------------------------------+
| | ExMin        | | Minimum of :math:`E_x`                                                  |
| | ExMinCell    | |  ... and its location (cell index)                                      |
| | ExMax        | | Maximum of :math:`E_x`                                                  |
| | ExMaxCell    | |  ... and its location (cell index)                                      |
| |              | | ... same for fields Ey Ez Bx_m By_m Bz_m Jx Jy Jz Rho                   |
| | PoyXmin      | | Accumulated Poynting flux through xmin boundary                         |
| | PoyXminInst  | | Current Poynting flux through xmin boundary                             |
| |              | |  ... same for other boundaries                                          |
+----------------+---------------------------------------------------------------------------+

Checkout the :doc:`post-processing <post-processing>` documentation as well.

----

.. _DiagFields:

*Fields* diagnostics
^^^^^^^^^^^^^^^^^^^^

:program:`Smilei` can collect various field data (electromagnetic fields, currents and density)
taken at the location of the PIC grid, both as instantaneous values and averaged values.
This is done by including a block ``DiagFields``::

  DiagFields(
      every = 10,
      time_average = 2,
      fields = ["Ex", "Ey", "Ez"]
  )

.. py:data:: every
  
  Number of timesteps between each output **or** a :ref:`time selection <TimeSelections>`.

.. py:data:: flush_every
  
  :default: 1
  
  Number of timesteps **or** a :ref:`time selection <TimeSelections>`.
  
  When `flush_every` coincides with `every`, the output
  file is actually written ("flushed" from the buffer). Flushing
  too often can *dramatically* slow down the simulation.


.. py:data:: time_average
  
  :default: ``1`` *(no averaging)*
  
  The number of timesteps for time-averaging.


.. py:data:: fields
  
  :default: ``[]`` *(all fields are written)*
  
  List of the field names that are saved. By default, they all are.


The full list of fields that are saved by this diagnostic:


.. rst-class:: nowrap

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
| | Rho          | |  Total density                                      |
| | Rho_abc      | |  Density of species "abc"                           |
+----------------+-------------------------------------------------------+


----

.. _DiagProbe:

*Probe* diagnostics
^^^^^^^^^^^^^^^^^^^

The fields from the previous section are taken at the PIC grid locations,
but it is also possible to obtain the fields at arbitrary locations.
These are called *probes*.

A probe interpolates the fields at either one point (0-D),
several points arranged in a line (1-D) or several points arranged in a mesh (2-D).

To add one probe diagnostic, include the block ``DiagProbe``::
  
  DiagProbe(
      every      = 10,
      pos        = [1., 1.],
      pos_first  = [1.,10.],
      pos_second = [10.,1.],
      number     = [100, 100],
      fields = ["Ex", "Ey", "Ez"]
  )

.. py:data:: every
  
  Number of timesteps between each output **or** a :ref:`time selection <TimeSelections>`.

.. py:data:: flush_every
  
  :default: 1
  
  Number of timesteps **or** a :ref:`time selection <TimeSelections>`.
  
  When `flush_every` coincides with `every`, the output
  file is actually written ("flushed" from the buffer). Flushing
  too often can *dramatically* slow down the simulation.


.. py:data:: pos
             pos_first
             pos_second
  
  :type: A list of floats, of length equal to the simulation dimensionality.
  
  | The coordinates of several points.
  | One point provided = a 0-D probe.
  | Two points provided = a 1-D probe.
  | Three points provided = a 2-D probe.

.. py:data:: number
  
  :type: A list of integers, one for each dimension of the probe.
  
  The number of points in each probe axis. Must not be defined for a 0-D probe.

.. py:data:: fields
  
  :default: ``[]`` (all fields)
  
  A list of fields among ``"Ex"``, ``"Ey"``, ``"Ez"``,
  ``"Bx"``, ``"By"``, ``"Bz"``, ``"Jx"``, ``"Jy"``, ``"Jz"`` and ``"Rho"``. Only these
  fields will be saved. 
  Note that it does NOT speed up calculation much, but it saves disk space.


**Examples of probe diagnostics**

* 0-D probe in 1-D simulation
  ::
    
    DiagProbe(
        every = 1,
        pos   = [1.2]
    )

* 1-D probe in 1-D simulation
  ::
    
    DiagProbe(
        every = 1,
        pos       = [1.2],
        pos_first = [5.6],
        number    = [100]
    )

* 1-D probe in 2-D simulation
  ::
    
    DiagProbe(
        every = 1,
        pos       = [1.2,  4.],
        pos_first = [5.6,  4.],
        number    = [100]
    )

* 2-D probe in 2-D simulation
  ::
    
    DiagProbe(
        every = 1,
        pos        = [0. ,   0.],
        pos_first  = [10. ,  0.],
        pos_second = [0.,    10.],
        number     = [100,   100]
    )


----

.. _DiagParticles:

*Particle* diagnostics
^^^^^^^^^^^^^^^^^^^^^^

A *particle diagnostic* collects data from the macro-particles and processes them during runtime.
It does not provide information on individual particles: instead, it produces
**averaged quantities** like the particle density, currents, etc.

The data is discretized inside a "grid" chosen by the user. This grid may be of any dimension.

Examples:

* 1-dimensional grid along the position :math:`x` (gives density variation along :math:`x`)
* 2-dimensional grid along positions :math:`x` and :math:`y` (gives density map)
* 1-dimensional grid along the velocity :math:`v_x` (gives the velocity distribution)
* 2-dimensional grid along position :math:`x` and momentum :math:`p_x` (gives the phase-space)
* 1-dimensional grid along the kinetic energy :math:`E_\mathrm{kin}` (gives the energy distribution)
* 3-dimensional grid along :math:`x`, :math:`y` and :math:`E_\mathrm{kin}` (gives the density map for several energies)
* 1-dimensional grid along the charge :math:`Z^\star` (gives the charge distribution)

Each dimension of the grid is called "axis".

You can add a particle diagnostic by including a block ``DiagParticles()`` in the namelist,
for instance::
  
  DiagParticles(
      output = "density",
      every = 5,
      time_average = 1,
      species = ["electrons1", "electrons2"],
      axes = [
          ["x", 0., 10, 100],
          ["ekin", 0.1, 100, 1000, "logscale", "edge_inclusive"]
      ]
  )

.. py:data:: output

  determines the data that is summed in each cell of the grid:
  
  * with ``"density"``, the weights are summed.
  * with ``"charge_density"``, the weights :math:`\times` charge are summed.
  * with ``"jx_density"``, the weights :math:`\times` charge :math:`\times\; v_x` are summed (same with :math:`y` and :math:`z`).
  * with ``"p_density"``, the weights :math:`\times\; p` are summed (same with :math:`p_x`, :math:`p_y` and :math:`p_z`).
  * with ``"ekin_density"``, the weights :math:`\times mc^2\; (\gamma-1)` are summed.
  * with ``"pressure_xx"``, the weights :math:`\times\; v_x p_x` are summed (same with yy, zz, xy, yz and xz).


.. py:data:: every
  
  The number of time-steps between each output, **or** a :ref:`time selection <TimeSelections>`.

.. py:data:: flush_every
  
  :default: 1
  
  Number of timesteps **or** a :ref:`time selection <TimeSelections>`.
  
  When `flush_every` coincides with `every`, the output
  file is actually written ("flushed" from the buffer). Flushing
  too often can *dramatically* slow down the simulation.


.. py:data:: time_average
  
  :default: 1
  
  The number of time-steps during which the data is averaged before output.


.. py:data:: species
  
  A list of the names of one or several species (see :py:data:`species_type`).


.. py:data:: axes
  
  A list of "axes" that define the grid.
  
  Syntax of one axis: ``[type, min, max, nsteps, "logscale", "edge_inclusive"]``
  
  * ``type`` is one of ``"x"``, ``"y"``, ``"z"``, ``"px"``, ``"py"``, ``"pz"``, ``"p"``,
    ``"gamma"``, ``"ekin"``, ``"vx"``, ``"vy"``, ``"vz"``, ``"v"`` or ``"charge"``.
  * The axis is discretized for ``type`` from ``min`` to ``max`` in ``nsteps`` bins.
  * The optional keyword ``logscale`` sets the axis scale to logarithmic instead of linear.
  * The optional keyword ``edge_inclusive`` includes the particles outside the range
    [``min``, ``max``] into the extrema bins.
  
  There may be as many axes as wanted in one ``DiagParticles( ... )`` block.

.. note::
  
  As an experimental capability, we created the "composite" axes ``type``.
  You may write the axis type as ``"ax+by+cz"``, where ``a``, ``b`` and ``c`` are numbers.
  This syntax does NOT accept characters other than numbers and the characters ``xyz+-``.
  For instance, it does not accept divisions ``/`` or whitespace.
  The resulting axis is along the vector of coordinates :math:`(a,b,c)`.
  For instance, in 2D, ``"x+2y"`` makes an axis oriented along the vector :math:`(1,2)`.


**Examples of particle diagnostics**

* Variation of the density of species ``electron1``
  from :math:`x=0` to 1, every 5 time-steps, without time-averaging
  ::
    
    DiagParticles(
    	output = "density",
    	every = 5,
    	time_average = 1,
    	species = ["electron1"],
    	axes = [ ["x",    0.,    1.,    30] ]
    )

* Density map from :math:`x=0` to 1, :math:`y=0` to 1
  ::
    
    DiagParticles(
    	output = "density",
    	every = 5,
    	time_average = 1,
    	species = ["electron1"],
    	axes = [ ["x",    0.,    1.,    30],
    	         ["y",    0.,    1.,    30] ]
    )

* Velocity distribution from :math:`v_x = -0.1` to :math:`0.1`
  ::
    
    DiagParticles(
    	output = "density",
    	every = 5,
    	time_average = 1,
    	species = ["electron1"],
    	axes = [ ["vx",   -0.1,    0.1,    100] ]
    )

* Phase space from :math:`x=0` to 1 and from :math:`px=-1` to 1
  ::
    
    DiagParticles(
    	output = "density",
    	every = 5,
    	time_average = 1,
    	species = ["electron1"],
    	axes = [ ["x",    0.,    1.,    30],
    	         ["px",   -1.,   1.,    100] ]
    )

* Energy distribution from 0.01 to 1 MeV in logarithmic scale.
  Note that the input units are :math:`m_ec^2 \sim 0.5` MeV
  ::
    
    DiagParticles(
    	output = "density",
    	every = 5,
    	time_average = 1,
    	species = ["electron1"],
    	axes = [ ["ekin",    0.02,    2.,   100, "logscale"] ]
    )

* :math:`x`-:math:`y` density maps for three bands of energy: :math:`[0,1]`, :math:`[1,2]`, :math:`[2,\infty]`.
  Note the use of ``edge_inclusive`` to reach energies up to :math:`\infty`
  ::
    
    DiagParticles(
    	output = "density",
    	every = 5,
    	time_average = 1,
    	species = ["electron1"],
    	axes = [ ["x",    0.,    1.,    30],
    	         ["y",    0.,    1.,    30],
    	         ["ekin", 0.,    6.,    3,  "edge_inclusive"] ]
    )

* Charge distribution from :math:`Z^\star =0` to 10
  ::
    
    DiagParticles(
    	output = "density",
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
It processes this data similarly to the :ref:`particle diagnostics <DiagParticles>`
as it makes a histogram of the macro-particle properties. The only difference is
that the histogram is made only by the particles that cross the surface.

You can add a screen by including a block ``DiagScreen()`` in the namelist,
for instance::
  
  DiagScreen(
      shape = "plane",
      point = [5., 10.],
      vector = [1., 0.],
      direction = "canceling",
      output = "density",
      species = ["electron"],
      axes = [["a", -10.*l0, 10.*l0, 40],
              ["px", 0., 3., 30]],
      every = 10
  )

.. py:data:: shape

   The shape of the screen surface: ``"plane"`` or ``"sphere"``.

.. py:data:: point

   :type: A list of floats ``[X]`` in 1D,  ``[X,Y]`` in 2D,  ``[X,Y,Z]`` in 3D 
   
   The coordinates of a point that defines the screen surface:
   a point of the ``"plane"`` or the center of the ``"sphere"``.

.. py:data:: vector

   :type: A list of floats ``[X]`` in 1D,  ``[X,Y]`` in 2D,  ``[X,Y,Z]`` in 3D 
   
   The coordinates of a vector that defines the screen surface:
   the normal to the ``"plane"`` or a radius of the ``"sphere"``.

.. py:data:: direction

   :default: ``"both"``
   
   Determines how particles are counted depending on which side of the screen they come from.
   
   * ``"both"`` to account for both sides.
   * ``"forward"`` for only the ones in the direction of the ``vector``.
   * ``"backward"`` for only the ones in the opposite direction.
   * ``"canceling"`` to count negatively the ones in the opposite direction.

.. py:data:: output

   Identical to the ``output`` of :ref:`particle diagnostics <DiagParticles>`.

.. py:data:: every
  
  The number of time-steps between each output, **or** a :ref:`time selection <TimeSelections>`.

.. py:data:: flush_every
  
  :default: 1
  
  Number of timesteps **or** a :ref:`time selection <TimeSelections>`.
  
  When `flush_every` coincides with `every`, the output
  file is actually written ("flushed" from the buffer). Flushing
  too often can *dramatically* slow down the simulation.

.. py:data:: species
  
  A list of the names of one or several species (see :py:data:`species_type`).

.. py:data:: axes
  
  A list of "axes" that define the grid of the histogram.
  It is identical to that of :ref:`particle diagnostics <DiagParticles>`, with the
  addition of four types of axes:
  ``"a"`` and ``"b"`` are the axes perpendicular to the ``vector``, when the screen
  shape is a ``"plane"``.
  ``"theta"`` and ``"phi"`` are the angles with respect to the ``vector``, when the screen
  shape is a ``"sphere"``.
  
----

.. _TimeSelections:

Time selections
^^^^^^^^^^^^^^^

Several components (mainly diagnostics) may require a selection of timesteps to
be chosen by the user. When one of these timesteps is reached, the diagnostics will
output data. A time selection is given through the parameter ``every`` and is a list
of several integers.

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

.. image:: _static/TimeSelections.png
  :width: 33em
  :align: center

..

.. admonition:: Tips
  
  * The syntax ``every = period`` is also accepted.
  * Any value set to ``0`` will be replaced by the default value.
  * Special case: ``every=0`` means no output.

----

.. _Checkpoints:

Checkpoints
^^^^^^^^^^^

The simulation can be *dumped* at given points (*checkpoints*) in order to be *restarted*
at that point.

A few things are important to know when you need dumps and restarts.

* Do not restart the simulation in the same directory as the previous one. Files will be 
  overwritten, and errors may occur. Create a new directory for your restarted simulation.
* Manage your memory: each process dumps one file, and the total can be significant.

::

  DumpRestart(
      restart_dir = "dump1",
      dump_step = 10000,
      dump_minutes = 240.,
      dump_deflate = 0,
      exit_after_dump = True,
      dump_file_sequence = 2,
  )

.. py:data:: restart_dir

  :default: None
  
  This tells :program:`Smilei` where to find dump files for restart.
  If not defined, it does not restart from a previous dump.
  
  **WARNING:** this path must either absolute or be relative to ``output_dir``
  
.. py:data:: dump_step

  :default: 0

  The number of timesteps between each dump of the full simulation.
  If ``0``, no dump is done.
  
.. py:data:: dump_minutes 

  :default: 0.

  The number of minutes between each dump of the full simulation (combines with ``dump_step``).
  If ``0.``, no dump is done.

.. py:data:: dump_deflate
  
  :red:`to do`

.. py:data:: exit_after_dump

  :default: ``True``

  If ``True``, the code stops after the dump.

.. py:data:: dump_file_sequence

  :default: 2
  
  This tells :program:`Smilei` to keep the last ``n`` dumps for a later restart.
  The default value, 2, saves one extra dump in case of a crash during the file dump.

.. py:data:: file_grouping

  :default: None
  
  The maximum number of checkpoint files that can be stored in one directory.
  New subdirectories are created according to the total number of files.
  This is useful on filesystem with a limited number of files per directory.

.. py:data:: restart_number

  :default: None
  
  If provided, the code will restart from that checkpoint, otherwise it uses the most recent one.

----

Variables defined by Smilei
^^^^^^^^^^^^^^^^^^^^^^^^^^^

:program:`Smilei` passes the following variables to the python interpreter for use in the
namelist. They should not be re-defined by the user!

.. py:data:: smilei_mpi_rank
    
  The MPI rank of the current process.

.. py:data:: smilei_mpi_size
    
  The total number of MPI processes.

.. py:data:: smilei_rand_max

  The largest random integer.


As an example of their use, this script randomizes both python's
and :program:`Smilei`'s random seeds.
::

    import random, math
    # reshuffle python random generator
    random.seed(random.random()*smilei_mpi_rank)
    # get 32bit pseudo random integer to be passed to smilei
    random_seed = random.randint(0,smilei_rand_max)
  

