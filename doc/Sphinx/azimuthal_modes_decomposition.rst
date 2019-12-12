
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

Indeed, the azimuthal decomposition concerns only the grid quantities (EM fields and current densities), but macro-particles evolve in a full three dimensional space.
Their positions and momenta are defined with 3D cartesian coordinates. 
At each iteration, they are evolved in the phase space as in a 3D simulation, using the 3D cartesian electromagnetic fields reconstructed from Eq. :eq:`AzimuthalDecomposition1`.
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

Defining diagnostics and initializing Profiles with s cylindrical geometry
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

If in doubt on how to initialize particles or a `Profile`, bear in mind how the quantities are defined in this geometry and the reference axes of the simulation (first Figure of this page).

Particles are defined in the 3D space, so if you want to initialize a `Species` with a numpy array you will still need to provide their coordinates
in the 3D cartesian space, but remembering that the :math:`y` and :math:`z` axes have their origins on the propagation axis (which is normally not the case in 2D and 3D cartesian simulations).
`Probes` diagnostics are like particles interpolating the reconstructed grid fields (including all the retained modes), so the same axes convention must be followed in defining their `origin` and `corners`.

Grid quantities instead are defined on the :math:`(x,r)` grid. Thus, `ExternalFields` and density/charge `Profiles` must be defined with functions of the :math:`(x,r)` coordinates.
Remember that `ExternalFields` are defined by mode. 





