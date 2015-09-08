Namelist syntax
---------------

In :program:`Smilei`, the input file (*namelist*) is written in the *python* language.
It is thus recommended to know the basics of *python*.

To create a namelist, we suggest you copy one existing file in the folder *benchmarks*.
All namelists have the extension *.py*.

----

General rules
^^^^^^^^^^^^^

* The namelist must define variables that :program:`Smilei` will understand.
  For instance, ``timestep = 0.01``.
  All *python* operations are valid. For instance: ``timestep = 40*0.0001``.

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
  

* You will need to use *lists*, which are series of things in *python*,
  defined between brackets ``[]`` and separated by commas.
  For example, ``mean_velocity = [0., 1.1, 3.]``.

* You are free to import any *python* package into the namelist.
  For instance, you may obtain :math:`\pi` using ``from math import pi``.

* **Important note**: :program:`Smilei` normalizes lengths and times to arbitrary values
  :math:`L_0` and :math:`T_0 = L_0/c`.
  All lengths and times defined in the namelist are given in units of
  :math:`L_0` and :math:`T_0`.

----

Stop and restart
^^^^^^^^^^^^^^^^
.. py:data:: dump_step 

  :default: 0

  The number of timesteps between each dump of the full simulation.
  If ``0``, no dump is done.
  
.. py:data:: dump_minutes 

  :default: 0.

  The number of minutes between each dump of the full simulation (combines with ``dump_step``).
  If ``0.``, no dump is done.

.. py:data:: exit_after_dump

  :default: ``True``

  If ``True``, the code stops after the dump.

.. py:data:: restart

  :default: ``False``

  If ``True``, :program:`Smilei` finds the last dump file and loads the corresponding simulation.
  If the dump file is not found, an error is raised.

.. py:data:: dump_file_sequence

  :default: 2
  
  This tells :program:`Smilei` to keep the last ``n`` dumps for a later restart 2 is the default option in case the code is stopped (or crashes) during a dump write leading to a unreadable dump file.
  
----

Spatial and temporal scales
^^^^^^^^^^^^^^^^^^^^^^^^^^^

.. py:data:: geometry
  
  The geometry of the simulation: ``"1d3v"`` or ``"2d3v"``.
  
  ``1d`` or ``2d`` correspond to the number of spatial dimensions.
  ``3v`` indicates the number of dimensions for velocities.


.. py:data:: interpolation_order
  
  :default: 2
  
  Interpolation order. To this day, only ``2`` is available.


.. py:data:: timestep
  
  Duration of one timestep in units of :math:`T_0`.


.. py:data:: sim_time
  
  Duration of the simulation in units of :math:`T_0`.


.. py:data:: cell_length
  
  A list of floats: dimensions of one cell in units of :math`L_0`.
  The number of elements of this list must be the same as the dimension of the simulation.


.. py:data:: sim_length
  
  A list of floats: dimensions of the simulations in units of :math:`L_0`.
  The number of elements of this list must be the same as the dimension of the simulation.


.. py:data:: clrw
  
  :default: 0.
  
  Cluster width.
  :red:`to do`


.. py:data:: wavelength_SI
  
  The value of the wavelength :math:`\lambda_0` in SI units
  (**only required if collisions or ionization are requested**).
  The wavelength is related to the normalization length according to :math:`2\pi L_0 = \lambda_0`.

.. py:data:: print_every
  
  Number of timesteps between each info output on screen. By default, 10 outputs per
  simulation.

----

Electromagnetic fields
^^^^^^^^^^^^^^^^^^^^^^

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


----

Moving window
^^^^^^^^^^^^^
.. py:data:: nspace_win_x

  :default: 0
  
  :red:`to do`


.. py:data:: t_move_win

  :default: 0.
  
  :red:`to do`


.. py:data:: vx_win

  :default: 0.
  
  :red:`to do`


----

Species
^^^^^^^
Each species has to be defined in a ``Species`` block, for instance::

  Species(
  	species_type = "electron",
  	initPosition_type = "regular",
  	initMomentum_type = "maxwell-juettner",
  	n_part_per_cell = 1000,
  	mass = 1.,
  	charge = 1.,
  	nb_density = 10.,
  	bc_part_type_west = "none",
  	bc_part_type_east = "none"
  )

All the possible variables inside this block are explained here:

.. py:data:: species_type
  
  :type: str
  
  The name you want to give to this species.


.. py:data:: initPosition_type
  
   The initialization of particle positions:
   
   * ``"regular"`` for regularly spaced
   * ``"random"`` for randomly distributed


.. py:data:: initMomentum_type
  
  The initialization of particle momenta:
  
  * ``"maxwell-juettner"`` for a relativistic maxwellian
  * ``"rectangular"`` for a rectangular distribution
  * ``"cold"`` for zero temperature
  
  The first 2 distributions depend on the parameter :py:data:`temperature` explained below.


.. py:data:: mass
  
  The mass of particles, in units of the electron mass :math:`m_e`.


.. py:data:: atomic_number

  The atomic number of the particles, required only if ionization is requested.
  :red:`todo`


.. py:data:: nb_density
             charge_density
  
  :type: float or *python* function (see section :ref:`profiles`)
  
  The absolute value of the number density or charge density (choose one only)
  of the particle distribution, in units of the critical density 
  :math:`n_c=\varepsilon_0 m_e/(e^2 T_0^2)`.


.. py:data:: charge
  
  :type: float or *python* function (see section :ref:`profiles`)
  
  The particle charge, in units of the electron charge $e$.


.. py:data:: mean_velocity
  
  :type: a list of 3 floats or *python* functions (see section :ref:`profiles`)
  
  The initial drift velocity of the particles, in units of the speed of light $c$.


.. py:data:: temperature
  
  :type: a list of 3 floats or *python* functions (see section :ref:`profiles`)
  
  The initial temperature of the particles, in units of :math:`m_ec^2`.


.. py:data:: n_part_per_cell
  
  :type: float or *python* function (see section :ref:`profiles`)
  
  The number of particles per cell.


.. py:data:: bc_part_type_west
             bc_part_type_east
             bc_part_type_south
             bc_part_type_north
  
  The boundary condition for particles: ``"none"`` means periodic.
  
  :red:`to do`


.. py:data:: time_frozen
  
  :default: 0.
  
  The time during which the particle positions are not updated, in units of :math:`T_0`.


.. py:data:: ionization_model
  
  :default: ``"none"``
  
  :red:`to do`


.. py:data:: radiating
  
  :default: ``False``
  
  :red:`to do`


.. py:data:: isTest
  
  :default: ``False``
  
  Flag for test particles. If ``True``, this species will contain only test particles
  which do not participate in the charge and currents.

.. py:data:: dump_every
  
  :default: 1
  
  Number of timesteps between each dump of test-particle information.
  Only active when ``isTest == True``


.. py:data:: c_part_max
  
  :red:`to do`


.. py:data:: dynamics_type
  
  :red:`to do`


----

External fields
^^^^^^^^^^^^^^^
External fields can be applied using a ``ExtField()`` block, for instance::

  ExtField(
      field = "Ex",
      profile = constant(0.01, xvacuum=0.1)
  )

All the possible variables inside this block are explained here:


.. py:data:: field
  
  The name of the field: ``"Ex"``, ``"Ey"``, ``"Ez"``, ``"Bx"``, ``"By"`` or ``"Bz"``.

.. py:data:: profile
  
  :type: float or *python* function (see section :ref:`profiles`)
  
  The initial spatial profile of the applied field.
  The units are the natural normalization units for fields,
  i.e., :math:`B_0=m_e /(eT_0)` for magnetic fields, and :math:`E_0=m_e c/(eT_0)` for electric fields.


----

.. _profiles:

Spatial profiles
^^^^^^^^^^^^^^^^

Several quantities require the input of a spatial profile:

* The average particle charge.
* The charge density or number density.
* The drift velocity.
* The temperature (only for distributions which require a temperature value).
* The number of particles per cell.
* Any external field.

There are many ways to define a profile.

1. Constant profiles

  * ``Species( ... , charge = -3., ... )`` defines a species with charge :math:`Z^\star=3`.
  
  * ``Species( ... , nb_density = 10., ... )`` defines a species with density :math:`10\,n_c`.
    You can choose ``nb_density`` (*number density*) or ``charge_density``
  
  * ``Species( ... , mean_velocity = [0.05, 0., 0.], ... )`` defines a species
    with drift velocity :math:`v_x = 0.05\,c` over the whole box.
  
  * ``Species(..., initMomentum_type="maxwell-juettner", temperature=[1e-5], ...)`` defines
    a species with a Maxwell-Jüttner distribution of temperature :math:`T = 10^{-5}\,m_ec^2` over the whole box.
    Note that the temperature may be anisotropic: ``temperature=[1e-5, 2e-5, 2e-5]``.
  
  * ``Species( ... , n_part_per_cell = 10., ... )`` defines a species with 10 particles per cell.
  
  * ``ExtField( field="Bx", profile=0.1 )`` defines a constant external field :math:`B_x = 0.1 B_0`.
  

2. *Python* profiles

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
  
..

  Once the function is created, you have to include it in the block you want:
  
  * ``Species( ... , charge = f, ... )``
  
  * ``Species( ... , nb_density = f, ... )``
  
  * ``Species( ... , mean_velocity = [f, f, f], ... )``
  
  * ``Species( ... , temperature = [f, f, f], ... )``
  
  * ``Species( ... , n_part_per_cell = f, ... )``
  
  * ``ExtField( ... , profile = f, ... )``


3. Built-in *python* functions

  :program:`Smilei` provides some *python* functions to help you build your profiles.
  
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
           xvacuum=0., xlength=None, phi=0., xnumber=1 )
  
    :param base: offset of the profile value
    :param amplitude: amplitude of the cosine
    :param xvacuum: empty length before starting the profile
    :param xlength: length of the profile (default is :py:data:`sim_length` :math:`-` ``xvacuum``)
    :param phi: phase offset`
    :param xnumber: number of periods within ``xlength``
  
  **Example**::
    
    Species( ... ,
        density = gaussian(10., xfwhm=0.3, xcenter=0.8),
    ... )

..
  
  **Illustration of the built-in profiles**
  
  .. image:: _static/pythonprofiles.png
  

----

Lasers
^^^^^^

:red:`to do`



----

.. _Collisions:

Collisions
^^^^^^^^^^

To have binary collisions in :program:`Smilei`,
add one or several ``Collisions()`` block in the namelist,
for instance::

  Collisions(
  	species1 = ["electrons1",  "electrons2"],
  	species2 = ["ions1"],
  	coulomb_log = 5.
  )


All the possible variables inside this block are explained here:

.. py:data:: species1
             species2
  
  :type: list of strings
  
  List of species names (see :py:data:`species_type`).
  
  The collisions will occur between
    1. all species under the list ``species1``
    2. and all species under the group ``species2``
  
  For instance, to have collisions between ``electrons1`` and ``ions1`` , use::
    
    species1 = ["electrons1"], species2 = ["ions1"]

..

  Other example, to collide all electrons with ions::
    
    species1 = ["electrons1", "electrons2"], species2 = ["ions"]

..

  **WARNING: this does not make** ``electrons1`` **collide with** ``electrons2``.
  
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


For more details about the collision scheme in :program:`Smilei`, see :doc:`collisions`


----

.. _DiagScalar:

*Scalars* diagnostics
^^^^^^^^^^^^^^^^^^^^^

:program:`Smilei` can collect various scalar data, such as total particle energy, total field energy, etc.
This is done by including the block ``DiagScalar()`` in the namelist, for instance::

  DiagScalar( every = 10 ,
              time_range = [0.1, 1.],
              vars = ["Utot", "Ukin", "Uelm"]
            )

All the possible variables inside this block are explained here:

.. py:data:: every
  
  Number of timesteps between each output.


.. py:data:: time_range
  
  :default: ``[]``
  
  | List of two values: minimum and maximum times that will be used.
  | Omit this argument to include all times.


.. py:data:: precision
  
  :default: 10
  
  Number of digits of the outputs.

.. py:data:: vars
  
  :default: ``[]``
  
  | List of scalars that will be actually output. Note that all scalars are computed anyways.
  | Omit this argument to include all scalars.


The full list of scalars that are saved by this diagnostic:

+----------------+---------------------------------------------------------------------------+
| **Global energies**                                                                        |
+----------------+---------------------------------------------------------------------------+
| | Utot         | | Total energy                                                            |
| | Ukin         | | Total kinetic energy (in the particles)                                 |
| | Uelm         | | Total EM energy (in the fields)                                         |
| | Uexp         | | Expected value (Initial energy :math:`-` lost :math:`+` gained)         |
| | Ubal         | | Energy balance (Utot :math:`-` Uexp)                                    |
| | Ubal_norm    | | Normalized energy balance (Ubal :math:`/` Utot)                         |
+----------------+---------------------------------------------------------------------------+
| **Energies lost/gained at boundaries due to moving window**                                |
+----------------+---------------------------------------------------------------------------+
| | Ukin_bnd     | | Kinetic energy exchanged at the boundaries during the timestep          |
| | Uelm_bnd     | | EM energy exchanged at boundaries during the timestep                   |
| | Ukin_out_mvw | | Kinetic energy lost during the timestep due to the moving window        |
| | Ukin_inj_mvw | | Kinetic energy injected during the timestep due to the moving window    |
| | Uelm_out_mvw | | EM energy lost during the timestep due to the moving window             |
| | Uelm_inj_mvw | | EM energy injected during the timestep due to the moving window         |
+----------------+---------------------------------------------------------------------------+
| **Energies lost/gained at boundaries  due to moving window**                               |
+----------------+---------------------------------------------------------------------------+
| | Ebal_norm    | | Ebalance :math:`/` Etot                                                 |
| | Ebalance     | | Current energy :math:`-` initial total energy                           |
| | Elost        | | Lost particle energy during last timestep                               |
| | Poynting     | | Accumulated Poyting flux through all boundaries                         |
+----------------+---------------------------------------------------------------------------+
| | Z_abc        | | Average charge of species "abc"                                         |
| | E_abc        | |  ... their kinetic energy                                               |
| | N_abc        | |  ... and number of particles                                            |
+----------------+---------------------------------------------------------------------------+
| | Ex_U         | | :math:`\int E_x^2 dV /2`                                                |
| |              | |  ... and similar for fields Ey, Ez, Bx_m, By_m and Bz_m                 |
+----------------+---------------------------------------------------------------------------+
| | ExMin        | | Minimum of :math:`E_x`                                                  |
| | ExMinCell    | |  ... and its location (cell index)                                      |
| | ExMax        | | Maximum of :math:`E_x`                                                  |
| | ExMaxCell    | |  ... and its location (cell index)                                      |
| |              | | ... same for fields Ey Ez Bx_m By_m Bz_m Jx Jy Jz Rho                   |
+----------------+---------------------------------------------------------------------------+
| | PoyEast      | | Accumulated Poynting flux through eastern boundary                      |
| | PoyEastInst  | | Current Poynting flux through eastern boundary                          |
| |              | |  ... same for boundaries West South North Bottom Top                    |
+----------------+---------------------------------------------------------------------------+

Checkout the :doc:`post-processing` documentation as well.

----

.. _DiagFields:

*Fields* diagnostics
^^^^^^^^^^^^^^^^^^^^

:program:`Smilei` can collect various field data (electromagnetic fields, currents and density)
taken at the location of the PIC grid, both as instantaneous values and averaged values.
This is done with the following instructions in the namelist:

.. py:data:: fieldDump_every
  
  The number of timesteps between each output of the instantaneous fields.

.. py:data:: avgfieldDump_every
  
  The number of timesteps between each output of the time-averaged fields.

.. py:data:: ntime_step_avg
  
  The number of timesteps for time-averaging.

.. py:data:: fieldsToDump
  
  :default: ``[]``
  
  List of field names that are saved. By default, they all are.


The full list of fields that are saved by this diagnostic:

+----------------+-------------------------------------------------------+
| | Bx_m         | |                                                     |
| | By_m         | | Components of the magnetic field                    |
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

To add one probe diagnostic, include the block ``DiagProbe()`` in the namelist.
There are several ways to do it:

**1. For only one point (zero-dimensional probe)**
  ::
    
    DiagProbe(
        every      = ... , # a number
        pos        = [x0, y0, z0]
    )
  
  * ``every`` is the number of timesteps between each output.
  * ``x0 [, y0 [, z0]]`` is the position of the point where to interpolate the fields.
  
  **Note**: ``y0`` (or ``z0``) should only be used in the case of a 2-D (or 3-D) simulation.


**2. For a series of points arranged in a line (one-dimensional probe)**
  ::
    
    DiagProbe(
        every      = ... , # a number
        pos        = [x0, y0, z0],
        pos_first  = [x1, y1, z1],
        number     = [n1]
    )
  
  * ``x0 [, y0 [, z0]]`` is the position of the starting point of the line.
  * ``x1 [, y1 [, z1]]`` is the position of the ending point of the line.
  * ``n1`` is the number of points along this line.

**3. For a series of points arranged in a mesh (two-dimensional probe)**
  ::
    
    DiagProbe(
        every      = ... , # a number
        pos        = [x0, y0, z0],
        pos_first  = [x1, y1, z1],
        pos_second = [x2, y2, z2],
        number     = [n1, n2]
    )
  
  In this case, the three points define three vertices of a paralellogram.


**Notes**

* There is an extra argument ``fields``, a list of fields among ``"Ex"``, ``"Ey"``, ``"Ez"``,
  ``"Bx"``, ``"By"``, ``"Bz"``, ``"Jx"``, ``"Jy"``, ``"Jz"`` and ``"Rho"``. Only these
  fields will be saved. Use, for example, ``fields=["Bz"]`` if you are only interested
  in :math:`B_z`. Note that it does NOT speed up calculation much, but it saves disk space.
* The dimension of the probe is decided only by the instruction ``number``:
  without it, the probe is 0-D, with ``number = [n1]``, the probe is 1-D,
  and with ``number =  [n1, n2]``, the probe is 2-D.
* You can have several probes in the input file.

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

The data may be collected from one or several particle species.

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

All the possible variables inside this block are explained here:

.. py:data:: output

  determines the data that is summed in each cell of the grid:
  
  * with ``"density"``, the weights are summed.
  * with ``"charge_density"``, the weights :math:`\times` charge are summed.
  * with ``"current_density_x"``, the weights :math:`\times` charge :math:`\times\; v_x` are summed (same with :math:`y` and :math:`z`).
  * with ``"p_density"``, the weights :math:`\times\; p` are summed (same with :math:`px`, :math:`py` and :math:`pz`)


.. py:data:: every
  
  The number of time-steps between each output.


.. py:data:: time_average
  
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


