
PIC algorithms
--------------

This document describes the theoretical basis with which :program:`Smilei` simulates
the behaviour of plasmas. Note that all quantities are expressed in terms of the
:doc:`reference units <units>`.

----

The Maxwell-Vlasov model
^^^^^^^^^^^^^^^^^^^^^^^^

The kinetic description of a collisionless plasma relies on the so-called Vlasov-Maxwell
system of equations. In this description, the different species of particles constituting
the plasma are described by their respective distribution functions
:math:`f_s(t,\mathbf{x},\mathbf{p})`, where :math:`s` denotes a given species consisting
of particles of charge :math:`q_s`, mass :math:`m_s`, and :math:`\mathbf{x}` and
:math:`\mathbf{p}` denote the position and momentum of a phase-space element.
The distribution :math:`f_s` satisfies Vlasov's equation:

.. math::
  :label: Vlasov

  \left(\partial_t + \frac{\mathbf{p}}{m_s \gamma} \cdot \nabla + \mathbf{F}_L \cdot \nabla_{\mathbf{p}} \right) f_s = 0\,,

where :math:`\gamma = \sqrt{1+\mathbf{p}^2/m_s^2}` is the (relativistic) Lorentz factor,

.. math::
  :label: LorentzForce

  \mathbf{F}_L = q_s\,(\mathbf{E} + \mathbf{v} \times \mathbf{B})

is the Lorentz force acting on the particles.
This force follows from the existence, in the plasma, of collective (macroscopic)
electric [:math:`\mathbf{E}(t,\mathbf{x})`] and magnetic [:math:`\mathbf{B}(t,\mathbf{x})`]
fields satisfying Maxwell's equations:

.. math::
  :label: Maxwell

  \begin{eqnarray}
  \nabla \cdot \mathbf{B} &=& 0 \,,\\
  \nabla \cdot \mathbf{E} &=& \rho \,,\\
  \nabla \times \mathbf{B} &=& \mathbf{J} + \partial_t \mathbf{E} \,,\\
  \nabla \times \mathbf{E} &=& -\partial_t \mathbf{B} \,.
  \end{eqnarray}

The Vlasov-Maxwell system of equations :eq:`Vlasov`--:eq:`Maxwell` describes the
self-consistent dynamics of the plasma which consistuents are subject to the Lorentz force,
and in turn modify the collective electric and magnetic fields through their charge and
current densities:

.. math::
  :label: rhoJ

  \begin{eqnarray}
  \rho(t,\mathbf{x}) &=& \sum_s q_s\int\!d^3\!p f_s(t,\mathbf{x},\mathbf{p})\,,\\
  \mathbf{J}(t,\mathbf{x}) &=& \sum_s q_s\int\! d^3\!p\,\mathbf{v} f_s(t,\mathbf{x},\mathbf{p})\,,
  \end{eqnarray}

where we have introduced the velocity :math:`\mathbf{v} = \mathbf{p}/(m_s\,\gamma)`.


----

.. _QuasiParticlesSection:

Quasi-particles
^^^^^^^^^^^^^^^

The *Particle-In-Cell* method owes its name to the discretization of the distribution
function :math:`f_s` as a sum of :math:`N_s` *quasi-particles* (also referred to as
*super-particles* or *macro-particles*):

.. math::
  :label: fs_discretized

  f_s(t,\mathbf{x},\mathbf{p}) =
    \sum_{p=1}^{N_s}\,w_p\,\,S\big(\mathbf{x}-\mathbf{x}_p(t)\big)\,\delta\big(\mathbf{p}-\mathbf{p}_p(t)\big)\,,

where :math:`w_p` is a quasi-particle *weight*, :math:`\mathbf{x}_p` is its position,
:math:`\mathbf{p}_p` is its momentum, :math:`S` is the shape-function of all quasi-particles,
and :math:`\delta` is the Dirac distribution.

In PIC codes, Vlasov's equation :eq:`Vlasov` is integrated along the continuous trajectories
of these quasi-particles, while Maxwell's equations :eq:`Maxwell` are solved on a
discrete spatial grid, the spaces between consecutive grid points being referred to as
*cells*. Injecting the discrete distribution function of
Eq. :eq:`fs_discretized` in Vlasov's equation :eq:`Vlasov`, multiplying the result by
:math:`\mathbf{p}` and integrating over all :math:`\mathbf{p}` and over the volume of
the quasi-particles, leads to the relativistic equations of motion of individual
quasi-particles:

.. math::

  \begin{eqnarray}
  \frac{d\mathbf{x}_p}{dt} &=& \frac{\mathbf{u}_p}{\gamma_p}\,\\
  \frac{d\mathbf{u}_p}{dt} &=& r_s \, \left( \mathbf{E}_p + \frac{\mathbf{u}_p}{\gamma_p} \times \mathbf{B}_p \right),
  \end{eqnarray}

where :math:`r_s = q_s/m_s` is the charge-over-mass ratio (for species :math:`s`),
:math:`\mathbf{u}_p = \mathbf{p}_p/m_s` is the reduced momentum and
:math:`\gamma_p=\sqrt{1+\mathbf{u}_p^2}` is the Lorentz factor.


----

Time and space discretization
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Maxwell's equations are solved here using
the `Finite Difference Time Domain (FDTD) approach <https://doi.org/10.1016/B978-012170960-0/50046-3>`_
as well as `refined methods based on this algorithm <https://doi.org/10.1140/epjd/e2014-50162-y>`_.
In these methods, the electromagnetic
fields are discretized onto a staggered grid, the so-called Yee-grid that allows for
spatial-centering of the discretized curl operators in Maxwell's equations.
The followingfigure summarizes at which points of the Yee-grid are defined the
electromagnetic fields as well as charge and density currents.

.. image:: _static/figYee.png
   :width: 13cm

Similarly, the time-centering
of the time-derivative in Maxwell's equations is ensured by considering the electric fields
as defined at integer time-steps :math:`(n)` and magnetic fields at half-integer
time-steps :math:`(n+\tfrac{1}{2})`. Time-centering of the magnetic fields is however
necessary for diagnostic purposes, and most importantly when computing the Lorentz force
acting on the quasi-particles.



A *leap-frog* scheme is used to advance the particles in time, so that the particle positions
and velocities are defined at integer :math:`(n)` and half-integer :math:`(n-\tfrac{1}{2})`
time-steps, respectively.

----

Initialization of the simulation
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

The initialization of a PIC simulation is a three-step process consisting in

#. loading particles,
#. computing the initial total charge and current densities on the grid,
#. computing the initial electric and magnetic field at the grid points.

In :program:`Smilei`, all three steps can be done either as a restart of a previous simulation
(in which case the particles, charge and current densities and electromagnetic fields are
directly copied from a file generated at the end of a previous simulation), or from a
user-defined input file. In that case, the user defines the initial conditions of the
particle, charge and current densities as well as the initial electromagnetic fields
over the whole simulation domain.

In particular, the number density :math:`n_s(\mathbf{x})`, mean velocity
:math:`\mathbf{v}_s(\mathbf{x})` and temperature :math:`T_s(\mathbf{x})` of all species
:math:`s` in a given cell (located at position :math:`\mathbf{x}`) at time :math:`t=0`
have to be prescribed. The particle loading then consists in creating, in each cell,
:math:`N_s` particles with positions :math:`\mathbf{x}_p` (either randomly chosen or
regularly spaced) such that particles are uniformly distributed within the cell,
and momentum :math:`\mathbf{p}_p` randomly chosen such that the particle distribution
follows a Maxwell-Jüttner distribution with mean-velocity :math:`\mathbf{v}_s(\mathbf{x})`
and temperature :math:`T_s(\mathbf{x})`.

In :program:`Smilei`, a weight is assigned to each particle depending on the density associated
to the cell it originates from:

.. math::

  w_p = \frac{n_s\big(\mathbf{x}_p(t=0)\big)}{N_s}\,.

This variable weighting is particularly beneficial when considering initially
highly-inhomogeneous density distributions.

Once all particles in the simulation domain have been created, the total charge and
current densities :math:`\rho(t=0,\mathbf{x})` and :math:`\mathbf{J}(t=0,\mathbf{x})`
are computed on the grid using a simple projection technique:

.. math::

  \rho(t=0,\mathbf{x}) = \sum_s\,q_s\,\sum_p\,w_p\,S\big(\mathbf{x}-\mathbf{x}_p(t=0)\big)\,.

Then, the initial electric fields are computed from :math:`\rho(t=0,\mathbf{x})`
by solving Poisson's equation. In :program:`Smilei`, this is done using the conjugate gradient
method. This iterative method is particularly interesting
as it is easily implemented on massively parallel computers and requires mainly
local information exchange between adjacent processes.

External (divergence-free) electric and/or magnetic fields can then be added to the
resulting electrostatic fields, provided they fullfill Maxwell's equations :eq:`Maxwell`,
and in particular Gauss' and Poisson's.

Note that a relativistic plasma needs :doc:`special treatment <relativistic_fields_initialization>`.

----

The PIC loop
^^^^^^^^^^^^

At the end of the initialization stage [time-step :math:`(n=0)`], all quasi-particles
in the simulation have been loaded and the electromagnetic fields have been computed
over the whole simulation grid. The PIC loop is then started over :math:`N` time-steps
each consisting in

#. interpolating the electromagnetic fields at the particle positions,
#. computing the new particle velocities and positions,
#. projecting the new charge and current densities on the grid,
#. computing the new electromagnetic fields on the grid.

In this section, we describe these four steps which advance the time from
time-step :math:`(n)` to time-step :math:`(n+1)`.


Field interpolation
"""""""""""""""""""

At the beginning of time-step :math:`(n)`, the particles velocity and position are known
at time-step :math:`(n-\tfrac{1}{2})` and :math:`(n)`, respectively. For each particle
:math:`p`, the electromagnetic fields [at time-step :math:`(n)`] are computed at the
particle position using a simple interpolation technique:

.. math::

  \begin{eqnarray}
  \mathbf{E}_p^{(n)} = V_c^{-1} \int d\mathbf{x}\, S\left(\mathbf{x}-\mathbf{x}_p^{(n)}\right) \mathbf{E}^{(n)}(\mathbf{x})\,,\\
  \mathbf{B}_p^{(n)} = V_c^{-1} \int d\mathbf{x}\, S\left(\mathbf{x}-\mathbf{x}_p^{(n)}\right) \mathbf{B}^{(n)}(\mathbf{x})\,,
  \end{eqnarray}

where we have used the time-centered magnetic fields
:math:`\mathbf{B}^{(n)}=\tfrac{1}{2}[\mathbf{B}^{(n+1/2) } + \mathbf{B}^{(n-1/2)}]`,
and :math:`V_c` denotes the volume of a cell.


Particle push
"""""""""""""

Knowing, for each quasi-particle, the electromagnetic fields at its position, the new
particle momentum and position are computed using a (second order) leap-frog integrator.

In :program:`Smilei`, different schemes have been implemented:
the well-known `Boris pusher <https://archive.org/stream/DTIC_ADA023511#page/n7/mode/2up>`_
both in the classical and relativistic form,
the `pusher developed by J.-L. Vay <https://doi.org/10.1063/1.2837054>`_,
and the `pusher of Higuera and Cary <https://arxiv.org/abs/1701.05605>`_.

All schemes compute the new particle momentum and position according to

.. math::

  \mathbf{u}_p^{n+\tfrac{1}{2}}=\mathbf{v}_p^{n-\tfrac{1}{2}} + r_s \Delta t \, \left[ E_p^{(n)} + \frac{\mathbf{v}_p^{(n+\tfrac{1}{2})}+\mathbf{v}_p^{(n-\tfrac{1}{2})}}{2} \times B_p^{(n)}\right],

.. math::

  \mathbf{x}_p^{n+1}=\mathbf{x}_p^{n} + \Delta t \, \frac{\mathbf{u}_p^{n+\tfrac{1}{2}}}{\gamma_p},

where :math:`\Delta t` denotes the duration of a time-step.


Current deposition
""""""""""""""""""

Charge deposition (i.e. charge and current density projection onto the grid) is then
performed using the charge-conserving algorithm
`proposed by Esirkepov <https://doi.org/10.1016/S0010-4655(00)00228-9>`_.
The current densities along the dimensions of the grid
(i.e., the :math:`x`-direction for 1D3V simulations,
both :math:`x`- and :math:`y`-directions for 2D3V simulations,
and all three :math:`x`-, :math:`y`- and :math:`z`-directions for 3D3V simulations)
are computed from the charge flux through the cell borders
(hence ensuring charge conservation) while the current densities along the other
dimensions are performed using a simple projection.

To illustrate this point, we take the example of current deposition in a 2D3V simulation.
The current densities in the :math:`x`- and :math:`y`-directions associated to a particle
with charge :math:`q` are computed as:

.. math::

  \begin{eqnarray}
  (J_x)_{i+\tfrac{1}{2},j}^{(n+\tfrac{1}{2})} = (J_x)_{i-\tfrac{1}{2},j}^{(n+\tfrac{1}{2})} + q\,w_p\,\frac{\Delta x}{\Delta t}\,(W_x)_{i+\tfrac{1}{2},j}^{(n+\tfrac{1}{2})}\,\\
  (J_y)_{i,j+\tfrac{1}{2}}^{(n+\tfrac{1}{2})} = (J_y)_{i,j-\tfrac{1}{2}}^{(n+\tfrac{1}{2})} + q\,w_p\,\frac{\Delta y}{\Delta t}\,(W_y)_{j,i+\tfrac{1}{2}}^{(n+\tfrac{1}{2})}\,
  \end{eqnarray}

where :math:`(W_x)^{(n+\tfrac{1}{2})}` and :math:`(W_y)^{(n+\tfrac{1}{2})}` are computed
from the particle current and former positions :math:`x_p^{(n+1)}` and :math:`x_p^{(n)}`,
respectively, using the method developed by Esirkepov.
The particle current in the :math:`z`-direction (not a dimension of the grid) is,
in this geometry, computed using a simple projection:

.. math::

  (J_z)_{i,j} = q w_r \mathbf{v}_p\,S(\mathbf{x}_{i,j}-\mathbf{x}_p)\,.


In all cases, the charge density deposited by the particle is obtained using the simple
projection:

.. math::

  (\rho)_{i,j}^{(n+1)} = q\,w_p\,S(\mathbf{x}_{i,j}-\mathbf{x}_p^{(n+1)})\,.

The total charge and current densities henceforth gather the contributions of all
quasi-particles of all species. It is worth noting that, within a charge-conserving
framework, charge densities are only projected on the grid for diagnostics purposes
(as we will see in the next paragraph, it is not used to advance the electromagnetic fields).


Maxwell solvers
"""""""""""""""

Now that the currents are known at time-step :math:`n+\tfrac{1}{2}`, the electromagnetic
fields can be advanced solving Maxwell's equations :eq:`Maxwell`.

First, Maxwell-Ampère is solved, giving the advanced electric fields

.. math::

  \mathbf{E}^{(n+1)} = \mathbf{E}^{(n)} + \Delta t\, \left[\left(\nabla \times \mathbf{B}\right)^{(n+\tfrac{1}{2})} - \mathbf{J}^{(n+\tfrac{1}{2})} \right]\,.

Then, Maxwell-Faraday is computed, leading to the advanced magnetic fields

.. math::

  \mathbf{B}^{(n+\tfrac{3}{2})} = \mathbf{B}^{(n+\tfrac{1}{2})} - \Delta t\, \left(\nabla \times \mathbf{E}\right)^{(n+1)}\,.

The discretization of the curl-operator is not detailed here.

It is worth
noting that computing the two previous equations is sufficient to get a complete description
of the new electromagnetic fields. Indeed, it can be shown that this conserves a
divergence-free magnetic field if Gauss' equation is satisfied at time :math:`t=0`.
Similarly, Poisson's equation is verified as long as it is satisfied
at time :math:`t=0`, if the charge deposition algorithm fulfills the charge conservation
equation:

.. math::

  \partial_t \rho + \nabla \cdot \mathbf{J} = 0

(this motivated the use of Esirkepov's projection scheme discussed in the previous paragraph).


----

Boundary conditions
^^^^^^^^^^^^^^^^^^^


After new quasi-particle positions and velocities have been computed, boundary conditions (BCs)
are applied to each quasi-particle that may be located in a ghost cell,
i.e. outside of the 'real' grid.
Quasi-particle species may have a different BC for each boundary of the simulation box:
the quasi-particles can either loop around the box (periodic),
be stopped (momentum set to zero),
suppressed (removed from memory),
reflected (momentum and position follow specular reflection rules)
or thermalized.
In the latter case, the quasi-particle is set back inside the simulation box,
and its new momentum is randomly sampled in a Maxwellian distribution
with a given temperature and drift velocity, both specified by the user.

BCs are applied to the electromagnetic fields after Maxwell's equations have been solved.
Each boundary of the simulation box can feature a different BC.
First, injecting/absorbing BCs inspired from the Silver-Müller BC
are able to inject an electromagnetic wave (e.g. a laser) and/or
to absorb outgoing electromagnetic waves.
In contrast, the reflective electromagnetic BC will reflect any outgoing
electromagnetic wave reaching the simulation boundary.
Lastly, periodic BCs correspond to applying the fields from the opposite boundary.


----

.. _multipassBinomialFilter:

Multi-pass binomial filtering of the current densities
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

A multi-pass binomial filter on the current densities is available in :program:`Smilei`,
which implementation follows that `presented by Vay et al. (2011) <https://www.sciencedirect.com/science/article/pii/S0021999111002270?via%3Dihub>`_.
Each pass consists in a 3-points spatial averaging (in all spatial dimensions) of the current densities, 
so that the filtered current density (here defined at location i on a one-dimensional grid) is recomputed as:

.. math::

    J_{f,i} = \frac{1}{2}\,J_i + \frac{J_{i+1}+J_{i-1}}{4}.


Current filtering, if required by the user, is applied before solving
Maxwell’s equation, and the number of passes is an :ref:`input parameter <CurrentFilter>`
defined by the user.



----

.. _EfieldFilter:

Friedman filter on the electric field
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

A method for temporal filtering of the electric field is also available in :program:`Smilei`.
It is the so-called Friedman filter detailed in `Greenwood et al. (2004) <https://www.sciencedirect.com/science/article/pii/S0021999104002608?via%3Dihub>`_.
This method consists in computing the filtered electric field at time-step :math:`n`:

.. math::

    {\bf E}_f^{(n)} = \left(1+\frac{\theta}{2}\right) {\bf E}^{(n)} - \theta \left(1-\frac{\theta}{2}\right) {\bf E}^{(n-1)} + \frac{1}{2} \theta \big(1-\theta\big)^2 \bar{\bf E}^{(n-2)},

where:

.. math::

    \bar{\bf E}^{(n-2)} = {\bf E}^{(n-2)} + \theta \bar{\bf E}^{(n-3)},

and :math:`\theta \in [0,1[` is an :ref:`input parameter <FieldFilter>` defined by the user.
Note that the filtered field :math:`E_f` is not used to push particles, but is used when solving the Maxwell-Faraday equation.
Also note that, as underlined in `Greenwood et al. (2004) <https://www.sciencedirect.com/science/article/pii/S0021999104002608?via%3Dihub>`_,
using this particular filter modifies the CFL condition of the Maxwell solver.
A simple trick to ensure that this condition is still verified is to use (for :math:`\Delta x = \Delta y = \Delta z`) the 
magic time-step :math:`\Delta t = \Delta x/2` whenever the Friedman filter is employed.

Both filters on the :ref:`currents <multipassBinomialFilter>` and :ref:`electric fields <EfieldFilter>` can be used together or
separately. They can be used, e.g., to mitigate the numerical Cherenkov instability that plagues PIC simulations dealing with 
relativistically drifting flows. 
An exemple of their use to mitigate this effect is highlighted in the work by `Plotnikov et al. (2017) <https://arxiv.org/abs/1712.02883>`_.


----

.. _AMcylindrical:

Azimuthal modes decomposition
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

:program:`Smilei` can run simulations in cyclindrical geometry with azimuthal modes decomposition as described in `this article <https://www.sciencedirect.com/science/article/pii/S0021999108005950?via%3Dihub>`_.
In this geometry, fields are expressed in cylindrical frame :math:`(e_x,e_r,e_\theta)` and can be written as a Fourier series:

.. math::

    F\left(x,r,\theta\right) = \textrm{Re}\left[\sum_{m=0}^{+\infty}\tilde{F}^{m}\left(x,r\right)\exp{\left(-im\theta\right)}\right],

where :math:`m` is the azimuthal mode, :math:`F` any field component and :math:`\tilde{F}^{m}` the :math:`m^{th}` Fourier mode of :math:`F`.
:math:`\tilde{F}^{m}` are given as follows:

.. math::

    \tilde{F}^{m} = \frac{1}{\pi}\int_0^{2\pi}F\left(x,r,\theta\right)\exp{\left(-im\theta\right)}d\theta,

for all :math:`m>0`. Mode 0 is given by:

.. math::

    \tilde{F}^{0} = \frac{1}{2\pi}\int_0^{2\pi}F\left(x,r,\theta\right)d\theta.

In the azimuthal modes decomposition simulations, only the :math:`\tilde{F}^{m}` are computed and stored.
Azimuthal modes from 0 to a maximum mode number defined by the user are considered.
These modes evolve independently of each other according to the linearity of Maxwell's equations.

From this collection of modes, a full three-dimensional reconstruction of the EM fields is possible using the equation above.
This can be done by using the :program:`Smilei` post processing tool :program:`Happi`. 
Note that each mode :math:`\tilde{F}^{m}` is a function of :math:`x`, the longitudinal coordinate and :math:`r`, the radial coordinate.
Therefore, each of them is only two dimensional.

Finally, note that this decomposition concerns only the grid quantities ( EM fields and current densities) and be aware that particles evolve in a full three dimensional space.


