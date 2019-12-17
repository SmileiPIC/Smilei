
Azimuthal modes decomposition
------------------------------------------

:program:`Smilei` can run simulations in cyclindrical geometry with azimuthal modes decomposition as described in `this article <https://www.sciencedirect.com/science/article/pii/S0021999108005950?via%3Dihub>`_.
Here we briefly review the method to give the necessary information to properly run a simulation in this geometry.

In many physical situations of interest, a propagation axis can be defined (e.g. referred to a laser pulse, a relativistic particle beam, ...). 
In :program:`Smilei` this axis is the `x` axis.

For a scalar field, e.g. a  density, perfect cylindrical symmetry around this propagation axis is defined as the invariance along the azimuthal angle :math:`\theta`, 
defined as in the following figure.

.. figure:: _static/Coordinate_Reference_AMcylindrical.png
  :width: 13cm
   
  Black axes: 3D cartesian reference for the particles coordinates and momenta. In blue the definition of the radial distance `r` and the angle :math:`\theta`


An azimuthal mode, or azimuthal harmonic, is defined as a complex exponential :math:`exp(-im\theta)`. 
Any scalar field can be decomposed into a basis of azimuthal modes.
The coordinates in this basis are Fourier coefficients and depend on the longitudinal and radial coordinates :math:`(x,r)`.
This azimuthal Fourier decomposition of a scalar field :math:`F(x,r,\theta)` would thus look like:

.. math::
    :label: AzimuthalDecomposition1

    F\left(x,r,\theta\right) = \textrm{Re}\left[\sum_{m=0}^{+\infty}\tilde{F}^{m}\left(x,r\right)\exp{\left(-im\theta\right)}\right],

where :math:`m` is the order of the azimuthal mode, :math:`\tilde{F}^{m}` the :math:`m^{th}` Fourier azimuthal mode of :math:`F`.

The perfectly cylindrical part of :math:`F`, constant along :math:`\theta` is given by the mode :math:`m=0`. 
The following figure shows the scalar fields defined through the real part of the azimuthal modes 
for :math:`m=0,1,2,3`. 


.. figure:: _static/AM_modes.png
  :width: 15cm
   
  Real part of the first four pure azimuthal modes :math:`exp(-im\theta)` on the `yz` plane.  

Expanding the series in Eq. :eq:`AzimuthalDecomposition1` yields:

.. math::
    :label: AzimuthalDecomposition2

    F\left(x,r,\theta\right) = \tilde{F}^{0}_{real} + \tilde{F}^{1}_{real}cos(\theta) + \tilde{F}^{1}_{imag}sin(\theta) + \tilde{F}^{2}_{real}cos(2\theta) + \tilde{F}^{2}_{imag}sin(2\theta) + ...


The complex coefficients :math:`\tilde{F}^{m}` of the azimuthal decomposition are given as follows for all modes of order :math:`m>0`:

.. math::

    \tilde{F}^{m>0} = \frac{1}{\pi}\int_0^{2\pi}F\left(x,r,\theta\right)\exp{\left(-im\theta\right)}d\theta,

The coefficient for the mode :math:`m=0` is given by:

.. math::

    \tilde{F}^{0} = \frac{1}{2\pi}\int_0^{2\pi}F\left(x,r,\theta\right)d\theta.

Note that also vector fields can be decomposed in azimuthal modes, through a decomposition of each of their components
along the cylindrical directions :math:`(e_x,e_r,e_\theta)`. 
For example, the transverse field :math:`\mathbf{E}_\perp` of a laser pulse polarized in the :math:`y` direction with cylindrically symmetric envelope
can be written as

.. math::

    \mathbf{E}_\perp(x,r,\theta, t) = E_y(x,r,\theta, t) e_y = E_r (x,r,\theta, t) e_r + E_{\theta}(x,r,\theta, t) e_{\theta} = E_y(x,r,t) [cos(\theta) e_r - sin(\theta) e_{\theta}].

Thus, referring to Eq :eq:`AzimuthalDecomposition2`, each of the cylindrical components of the mentioned laser at a given instant would be composed of a pure azimuthal mode of order :math:`m=1`, 
multiplied by its Fourier coefficient :math:`\tilde{E}^1(x,r,t)`:

.. math::

    \tilde{E}^1_r (x,r,\theta) = E_y(x,r,t),\\

    \tilde{E}^1_{\theta} (x,r,\theta) = -iE_y(x,r,t).

Similarly, an elliptically (or cilindrically) polarized laser would be given by an azimuthal decomposition of their cylindrical components,
with only the mode :math:`m=1`. Indeed, a laser with elliptical polarization can be seen as the linear superposition of two linearly polarized lasers,
with different phases and amplitudes. The difference in phase would be equivalent to the multiplication of the Fourier coefficient by a complex exponential.

Physical phenomena close to cylindrical symmetry, where the use of simulations with this technique is most suited, can in principle be characterised only by the presence of
the low order azimuthal modes, since the Fourier coefficients of the higher order modes (representing stronger cylindrical asymmetry) are zero or negligible.

For example, in a basic Laser Wakefield Acceleration setup, a laser pulse with cylindrically symmetric envelope could be described only by the mode :math:`m=1` and the cylindrically symmetric wave
in its wake by the mode :math:`m=0`. This because the shape of the wake wave is mainly determined by the ponderomotive force of the laser, which depends on its cylindrically symmetric envelope. 
Thus, a simulation of this phenomenon would in principle need only two azimuthal modes. In the namelist of the corresponding simulation with azimuthal modes decomposition (`geometry=AMcylindrical`), 
the user would then choose `number_of_AM=2` in this case.  

In the azimuthal modes decomposition simulations, only the :math:`\tilde{F}^{m}` of Eq. :eq:`AzimuthalDecomposition1`, for each scalar field and for all the components of the vector fields, 
are computed and stored. Each of them is a complex field defined in the :math:`(x,r)` space.
In other words, for all the physical grid fields only the azimuthal modes from 0 to `(number_of_AM-1)` are considered. 

In vacuum, the azimuthal modes of the cylindrical components of the electromagnetic fields would evolve independently. 
Due to the linearity of Maxwell's Equations, we can write and solve them separately for each mode.
The resulting equations describing the mode :math:`m` evolution in presence of current densities are:

.. math::
    :label: MaxwellEqsAzimuthalModes

    \partial_t \tilde{B}^m_{x} =-\frac{1}{r}\partial_r(r\tilde{E}^m_{\theta})-\frac{im}{r}\tilde{E}^m_r,\\
    \partial_t \tilde{B}^m_r = \frac{im}{r}\tilde{E}^m_x+\partial_x \tilde{E}^m_{\theta},\\
    \partial_t \tilde{B}^m_{\theta} =-\partial_x \tilde{E}^m_{r} + \partial_r \tilde{E}^m_{x},\\
    \partial_t \tilde{E}^m_{x} =\frac{1}{r}\partial_r(r\tilde{B}^m_{\theta})+\frac{im}{r}\tilde{B}^m_r-\tilde{J}^m_{x},\\
    \partial_t \tilde{E}^m_r = -\frac{im}{r}\tilde{B}^m_x-\partial_x \tilde{B}^m_{\theta}-\tilde{J}^m_{r},\\
    \partial_t \tilde{E}^m_{\theta} =\partial_x \tilde{B}^m_{r} - \partial_r \tilde{B}^m_{x}-\tilde{J}^m_{\theta}.

Thus even in presence of a plasma (i.e. non zero current densities), at each timestep these equations are solved independently. 
The coupling between the modes occurs when the electromagnetic fields (the superposition of their the modes obtained from Eq. :eq:`AzimuthalDecomposition1`) interact with the particles, 
which in turn create the sources for Eqs. :eq:`MaxwellEqsAzimuthalModes`, i.e. the azimuthal components :math:`\tilde{J}^m` of their current density.  

Indeed, the azimuthal decomposition concerns only the grid quantities (EM fields and current densities), which are thus defined on a 2D grid, but macro-particles evolve in a full three dimensional space.
Their positions and momenta are defined with 3D cartesian coordinates. This concept is illustrated in the following figure.

.. figure:: _static/AM_grid_particles.png
  :width: 10cm
   
  Definition of the grid quantities and of the particles positions in a simulation with azimuthal modes decomposition. Blue arrows: the `x` and `r` axes of the 2D grid (in red) where the electromagnetic fields are defined. The particles positions and momenta are defined in the 3D space.

At each iteration, the particles are evolved in the phase space as in a 3D simulation, using the 3D cartesian electromagnetic fields reconstructed from Eq. :eq:`AzimuthalDecomposition1`.
The angle :math:`\theta` of each particle is computed from its position to know the total electromagnetic field acting on it (reconstructed from Eq. :eq:`AzimuthalDecomposition1` ).
Then, again depending on their position :math:`(x,r,\theta)`, their azimuthal contribution to the current densities :math:`(J_x,J_r,J_{\theta})` are computed 
to evolve the electromagnetic fields at the next PIC iteration solving Eqs :eq:`MaxwellEqsAzimuthalModes`.

The same reconstruction can be done through the :program:`Smilei` post processing tool :program:`Happi`.  
Note that each mode :math:`\tilde{F}^{m}` is a function of :math:`x`, the longitudinal coordinate and :math:`r`, the radial coordinate.
Therefore, each of them is only two dimensional. Thus, the computational cost of simulations with azimuthal decompositions in principle scales approximately as 
`number_of_AM` simulations in 2D, but obtaining results with 3D accuracy if a suitable number of modes is used. 
Note although that, due to the cylindrical geometry, a higher number of particles than in a 2D  or 3D cartesian simulation could be necessary to obtain convergence of the results.
In this geometry, always check the convergence of your results trying to increase the number of macro-particles and of retained modes.
A rule of thumb is to use at least :math:`4\times number\_of\_AM` macro-particles along :math:`\theta`.


----

Defining diagnostics and initializing Profiles with a cylindrical geometry
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

If in doubt on how to initialize particles or a `Profile`, bear in mind how the quantities are defined in this geometry and the reference axes of the simulation (first Figure of this page).

Note also in the following figure the difference in the origin of the reference axes between a simulation with azimuthal modes decomposition and a 3D Cartesian simulation.
In the figure, the red ellipsoid can represent a laser pulse or a relativistic particle beam, whose propagation axis is the red dashed arrow. In the case of the simulation with azimuthal modes, this propagation axis is the `x` axis. In a 3D simulation for the same physical case,
this axis would be parallel to the `x` axis. If this is not the case, probably the azimuthal modes decomposition technique is not suited for the simulation, since too many azimuthal modes
would be necessary for an accurate representation.

.. figure:: _static/AMcylindrical_vs_cartesian.pdf
  :width: 22cm
   
  Comparison between a simulation with azimuthal modes decomposition and a 3D Cartesian simulation for the same physical case. The blue point is the origin of the reference axes. 
  The radial grid size (`grid_length[1]`) of the former is half the size of the `y-z` grid sizes (`grid_length[1]=grid_length[2]`) of the latter.

Particles are defined in the 3D space, so if you want to initialize a `Species` with a numpy array you will still need to provide their coordinates
in the 3D cartesian space.
`Probes` diagnostics with azimuthal modes are like particles interpolating the reconstructed grid fields (including all the retained modes), so the same axes convention of the previous figure must be followed in defining their `origin` and `corners`.

Grid quantities instead are defined on the :math:`(x,r)` grid. Thus, `ExternalFields` and density/charge `Profiles` must be defined with functions of the :math:`(x,r)` coordinates.
Remember that `ExternalFields` are defined by mode. 



----

Poisson's equation and relativistic Poisson's equation with azimuthal modes decomposition
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

In a simulation with azimuthal modes decomposition, given the linearity of the relativistic Poisson's equation described in :doc:`relativistic_fields_initialization`, the full equation
can be decomposed in azimuthal modes, with the correspondent mode component of the charge density :math:`-\tilde{\rho}^m` as source term.

The relativistic Poisson equation for the potential component :math:`\tilde{\Phi}^m` of the mode :math:`m` in this  geometry is thus:

.. math::
  :label: RelPoissonModes

  \left[ \frac{1}{\gamma^2_0}\partial^2_x\tilde{\Phi}^m+\frac{1}{r}\partial_r\left(r\partial_r\tilde{\Phi}^m\right)-\frac{m^2}{r^2}\tilde{\Phi}^m \right] = -\tilde{\rho}^m.

Solving each of these relativistic Poisson's equations allows to initialize the azimuthal components of the electromagnetic fields:

.. math::
  \begin{eqnarray} 
  \tilde{E}^m_x &=& -\frac{1}{\gamma_0^2}\partial_x \tilde{\Phi}^m,\\ 
  \tilde{E}^m_r &=& -\partial_r \tilde{\Phi}^m, \\ 
  \tilde{E}^m_{\theta} &=& \frac{im}{r} \tilde{\Phi}^m,\newline\\
  \tilde{\mathbf{B}}^m &=& \beta_0\mathbf{\hat{x}}\times\tilde{\mathbf{E}}^m.
  \end{eqnarray} 

The initialization of the electric field with the non relativistic Poisson's equation is performed similarly, and the underlying equations reduce simply
to the previous equations, but with :math:`\gamma_0 = 1` and :math:`\beta_0 = 0` (i.e. an immobile Species).


----

The envelope model in cylindrical coordinates
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

In :program:`Smilei` the :doc:`laser_envelope` described in [Terzani]_, [MassimoPPCF2019]_ for cartesian geometries has been implemented also in cylindrical geometry,
as described in [Massimo2020]_.

The azimuthal decomposition technique is used in this case, but only the mode :math:`m=0` can be retained in the present implementation, 
i.e. the electromagnetic fields and the envelope fields will have perfect cylindrical symmetry with respect to the envelope propagation axis :math:`x`.

The main difference compared to the cartesian geometry lies in the envelope equation, Eq. :eq:`envelope_equation` of the page :doc:`laser_envelope`. 
While the Laplacian operator :math:`\nabla^2` is defined as
:math:`\partial_x^2`, :math:`\partial_x^2+\partial_y^2` and :math:`\partial_x^2+\partial_y^2+\partial_z^2` 
in 1D, 2D and 3D cartesian coordinates respectively, the envelope equation in `AMcylindrical` geometry of course uses the Laplacian in 
cylindrical coordinates. Additionally, due to the assumption of cylindrical symmetry, the derivatives with respect to the azimuthal angle are all zero by definition.
Thus, in this geometry the envelope equation solved in :program:`Smilei` is:

.. math::
  :label: envelope_equation

  \partial^2_x\tilde{A}+\frac{1}{r}\partial_r(r\partial_r\tilde{A})+2i\left(\partial_x \tilde{A} + \partial_t \tilde{A}\right)-\partial^2_t\tilde{A}=\chi \tilde{A}.

The electromagnetic fields evolve as described in :doc:`azimuthal_modes_decomposition` with only the mode :math:`m=0`, 
or equivalently neglecting all the derivatives along the azimuthal angle in Maxwell's Equations written in cylindrical coordinates.

As in a typical :program:`Smilei` simulation in cylindrical coordinates, the particles evolve in the 3D space, 
with their positions and momenta described in cartesian coordinates.

The envelope approximation coupled to the cylindrical symmetry assumption can greatly speed-up a simulation of a physical set-up where these assumptions are suited.
Compared to a 3D envelope simulation with the same number of particles, a cylindrical envelope simulation has a speed-up which scales linearly 
as the double of the transverse number of cells of the window. This speed-up can arrive to at least a factor 100 for lasers with transverse sizes of the order of tens of microns.
Compared to a 3D standard laser simulation with the same number of particles, 
the speed-up of a cylindrical envelope simulation can arrive to at least a factor 1000 for lasers of durations of the order of tens of femtoseconds. 
These comparisons assume the same longitudinal window size and the same transverse size for the simulated physical space.




