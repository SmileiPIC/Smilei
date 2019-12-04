
Azimuthal modes decomposition
------------------------------------------

:program:`Smilei` can run simulations in cyclindrical geometry with azimuthal modes decomposition as described in `this article <https://www.sciencedirect.com/science/article/pii/S0021999108005950?via%3Dihub>`_.
Here we briefly review the method to give the necessary information to properly run a simulation in this geometry.

In many physical situations of interest, a propagation axis can be defined (e.g. referred to a laser pulse, a relativistic particle beam, ...). 
In :program:`Smilei` this axis is the `x` axis.

For a scalar field, e.g. a  density, perfect cylindrical symmetry referred to this propagation axis is defined as the invariance along the azimuthal angle :math:`\theta`, 
defined as in the following figure.

.. figure:: _static/Coordinate_Reference_AMcylindrical.png
  :width: 13cm
   
  Black axes: 3D cartesian reference for the particles coordinates and momenta. In blue the definition of the radial distance `r` and the angle :math:`\theta`

We can define an azimuthal mode, or azimuthal harmonic as a complex exponential :math:`exp(-im\theta)`. 
The perfect cylindrical symmetry corresponds to :math:`m=0`, a constant. For the sake of illustration, in the following figure the scalar fields defined through the real part of the azimuthal modes 
for :math:`m=0,1,2,3` are shown. 


.. figure:: _static/AM_modes.png
  :width: 15cm
   
  Real part of the first four pure azimuthal modes :math:`exp(-im\theta)` on the `yz` plane, with no radial dependence.  


If we decomposed a physical scalar field in azimuthal modes, the Fourier coefficients for the azimuthal modes would depend also
on the longitudinal and radial coordinates :math:`(x,r)`. The azimuthal Fourier decomposition of a scalar field :math:`F(x,r,\theta)` would thus look like:

.. math::
    :label: AzimuthalDecomposition

    F\left(x,r,\theta\right) = \textrm{Re}\left[\sum_{m=0}^{+\infty}\tilde{F}^{m}\left(x,r\right)\exp{\left(-im\theta\right)}\right],

where :math:`m` is the azimuthal mode :math:`\tilde{F}^{m}` the :math:`m^{th}` Fourier mode of :math:`F`.
The coefficients :math:`\tilde{F}^{m}` of the azimuthal decomposition are given as follows:

.. math::

    \tilde{F}^{m} = \frac{1}{\pi}\int_0^{2\pi}F\left(x,r,\theta\right)\exp{\left(-im\theta\right)}d\theta,

for all modes of order :math:`m>0`. The coefficient for the mode :math:`m=0` is given by:

.. math::

    \tilde{F}^{0} = \frac{1}{2\pi}\int_0^{2\pi}F\left(x,r,\theta\right)d\theta.

Note that also vector fields can be decomposed in azimuthal modes, through a decomposition of each of their components
along the cylindrical directions :math:`(e_x,e_r,e_\theta)`. 
For example, the transverse field :math:`\mathbf{E}_\perp` of a laser pulse polarized in the :math:`y` direction with cylindrically symmetric envelope
can be written as

.. math::

    \mathbf{E}_\perp(x,r,\theta, t) = E_y(x,r,\theta, t) e_y = E_r (x,r,\theta, t) e_r + E_{\theta}(x,r,\theta, t) e_{\theta} = E_y(x,r,t) [cos(\theta) e_r - sin(\theta) e_{\theta}].

Thus, referring to Eq :eq:`AzimuthalDecomposition`, each of the cylindrical components of the mentioned laser at a given instant would be composed of a pure azimuthal mode of order :math:`m=1`, 
multiplied by its Fourier coefficient :math:`\tilde{E}^1(x,r,t)`:

.. math::

    \tilde{E}^1_r (x,r,\theta) = E_y(x,r,t),\\

    \tilde{E}^1_{\theta} (x,r,\theta) = -iE_y(x,r,t).

Similarly, an elliptically (or cilindrically) polarized laser would be given by an azimuthal decomposition of their cylindrical components,
with only the mode :math:`m=1`. Indeed, a laser with elliptical polarization can be seen as the linear superposition of two linearly polarized lasers,
with different phases and amplitudes. The difference in phase would be equivalent to the multiplication of the Fourier coefficient by a complex exponential.

Physical phenomena with a high degree cylindrical symmetry, where the use of simulations with this technique is most suited, can in principle be characterised only by the presence of
the low order azimuthal modes, since the Fourier coefficients of the higher order modes (representing a high degree of cylindrical asymmetry) are zero or negligible.

For example, in a basic Laser Wakefield Acceleration set-up a laser pulse with cylindrically symmetric envelope could be described only by the mode :math:`m=1` and the cylindrically symmetric wave
in its wake by the mode :math:`m=0`. Thus, a simulation of this phenomenon would need only two azimuthal modes (`number_of_AM=2` in the input namelist). 

In the azimuthal modes decomposition simulations, only the :math:`\tilde{F}^{m}` of Eq. :eq:`AzimuthalDecomposition`, for each scalar field and for all the components of the vector fields, 
are computed and stored. Each of them is a complex field defined in the :math:`(x,r)` space.
In other words, for all the physical grid fields only the azimuthal modes from 0 to `(number_of_AM-1)`, the latter parameter defined by the user in the namelist, are considered. 

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

Even in presence of a plasma, at each timestep these equations are solved independently. 
The coupling between the modes occurs when the electromagnetic fields (the superposition of their the modes) interact with the particles, 
which in turn create the sources for Eqs. :eq:`MaxwellEqsAzimuthalModes`, i.e. the azimuthal components :math:`\tilde{J}^m` of their current density.  

Indeed, the azimuthal decomposition concerns only the grid quantities (EM fields and current densities), but particles evolve in a full three dimensional space.
Their positions and momenta are defined with 3D cartesian coordinates. 
At each iteration, they are evolved in the phase space as in a 3D simulation, using the 3D cartesian electromagnetic fields reconstructed from Eq. :eq:`AzimuthalDecomposition`.
Then, depending on their position :math:`(x,r,\theta)`, their azimuthal contribution to the current densities :math:`(J_x,J_r,J_{\theta})` are computed 
to evolve the electromagnetic fields at the next PIC iteration solving Eqs :eq:`MaxwellEqsAzimuthalModes`.

The same reconstruction can be done through the :program:`Smilei` post processing tool :program:`Happi`.  
Note that each mode :math:`\tilde{F}^{m}` is a function of :math:`x`, the longitudinal coordinate and :math:`r`, the radial coordinate.
Therefore, each of them is only two dimensional. The computational cost of simulations with azimuthal decompositions in principle scales then approximately as 
`number_of_AM` simulations in 2D, but obtaining results with 3D accuracy if a suitable number of modes is used. 



