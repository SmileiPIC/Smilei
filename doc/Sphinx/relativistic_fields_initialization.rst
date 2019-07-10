Field initialization for relativistic species
--------------------------------------------------------------------------------

As explained in :doc:`algorithms`, if a net charge is present at the beginning of the simulation, the initial electromagnetic fields are computed.
For static charge distributions, the solution of Poisson's equation will be necessary to find the initial electrostatic field. 
If the initial charge has a non-zero initial speed, in general the electric and magnetic field should be computed solving the full set of Maxwell's equations or equivalently the potentials equations.
In some physical setups of interest, one or more relativistic species are injected in a plasma. In these cases, the computation of the initial electromagnetic fields can be reduced to the solution of a modified version of Poisson's equation.


----

The relativistic Poisson's equation
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

From the continuity equation :math:`\partial_t \rho + \nabla \cdot \mathbf{J} = 0`
and Maxwell's equations, it can be shown that the quantities :math:`\nabla\cdot\mathbf{B}` and :math:`\nabla\cdot\mathbf{E}-\rho` do not change over time:

.. math::

  \begin{eqnarray}
  \partial_t \left( \nabla\cdot\mathbf{B} \right ) &=& 0, \\
  \partial_t \left( \nabla\cdot\mathbf{E}-\rho \right ) &=& \nabla\cdot\partial_t\mathbf{E}-\partial_t\rho = \nabla\cdot\left(\nabla\times\mathbf{B}-\mathbf{J}\right)-\partial_t\rho = - \left(\nabla\cdot\mathbf{J}+\partial_t \rho\right).
  \end{eqnarray}

Thus, if a simulation starts with :math:`\rho\neq0`, the electromagnetic fields must be properly initialized. 

In the case of a static charge distribution, i.e. :math:`\rho\neq0`, :math:`\mathbf{J}=0`, the initial electrostatic potential :math:`\Phi` can be computed solving Poisson's equation:

.. math::

  \nabla^2 \Phi = -\rho,

and then be integrated to find the initial electric field: :math:`\mathbf{E}=-\nabla\Phi`. The initial magnetic field :math:`\mathbf{B}` will be zero.

In general when the initial current :math:`\mathbf{J}` is not zero, the full set of fields equations should be solved to correctly initialize the electromagnetic fields. 

However, if a species is already relativistic when the simulation starts, e.g. a relativistic electron bunch, its initial electromagnetic fields can be computed through a simplified procedure, described in [Vay2008]_, [Londrillo2014]_, [Massimo2016]_ and [Marocchino2018]_. 

An important assumption of this calculation is that the species is highly relativistic, moving in the positive :math:`x` direction, with negligible momentum spread. Under this hypothesis, the transverse components of the species current density are neglected and the four-current quadrivector can be written as:

.. math::

  \left(\mathbf{J},\rho\right) = \left(\rho \beta_0 , 0, 0, \rho \right),

where :math:`\beta_0` is the initial mean velocity of the relativistic species. At least locally, the potentials :math:`\mathbf{A}`, :math:`\Phi` in the laboratory frame will be only function of :math:`x-\beta_0 t`, as they are propagating with the species at uniform relativistic velocity.

In the relativistic species rest frame :math:`S'`, the charge distribution is static and the electrostatic potential in that reference frame :math:`\Phi'` is related to the charge density in that reference frame :math:`\rho'` through Poisson's equation:

.. math::
  :label: Poisson

  \nabla'^2 \Phi' = -\rho',

where the Laplacian operator is computed in the reference frame :math:`S'`:

.. math::
  
  \nabla'^2=\partial^2_{x'}+\partial^2_{y'}+\partial^2_{z'}.

The vector potential in the species rest frame can be set to zero: :math:`\mathbf{A'}=0`. Through the above mentioned assumptions, it is possible to rewrite Eq. :eq:`Poisson` only in terms of laboratory frame quantities. 

Lorentz transformation of the four-vector :math:`\left(\mathbf{J},\rho \right)` yields :math:`\rho'=\rho/\gamma_0`, where :math:`\gamma_0=1/\sqrt{1-\beta^2_0}` is the average Lorentz factor of the relativistic species. 
Similarly, the potential :math:`\Phi'` can be rewritten in terms of the potential in the laboratory frame: :math:`\Phi'=\Phi/\gamma_0`. The Lorentz back-transformation of coordinates

.. math::
  
  x=\gamma_0(x'+\beta_0 t'),\quad  t = \gamma_0(t'+\beta_0 x'), \quad y=y', \quad z=z'

allows to transform the derivatives in Eq. :eq:`Poisson` as 

.. math::
  
  \partial_{x'}=\gamma_0\left(\partial_x+\beta_0\partial_t\right), \quad \partial_{y'}=\partial_y, \quad \partial_{z'}=\partial_z. 

The partial derivative along the :math:`x'` direction can be further simplified, through the hypothesis of temporary dependence of all quantities on :math:`x-\beta_0 t`, implying :math:`\partial_t=-\beta_0\partial_x`:

.. math::
  
  \partial_{x'}=\frac{1}{\gamma_0}\partial_x. 

Equation :eq:`Poisson` can thus be rewritten as 

.. math::
  :label: RelPoisson

  \left( \frac{1}{\gamma^2_0}\partial^2_x+\nabla_{\perp}^2\right) \Phi = -\rho,

here informally referred to as the relativistic Poisson's equation. In :program:`Smilei`, as for Eq. :eq:`Poisson`, the solution of the relativistic Poisson's equation is performed through the conjugate gradient method.

Once the potential :math:`\Phi` is found, we can compute all the components of the electromagnetic field, using again the relations :math:`\partial_t=-\beta_0\partial_x`, :math:`\Phi'=-\Phi/\gamma_0` and the Lorentz back-transformation of the vector potential :math:`\mathbf{A}`:

.. math::
  
  A_x = \gamma_0(A_x'+\beta_0 \Phi')=\gamma_0\beta_0 \Phi'=\beta_0\Phi,\quad A_y = A_y'=0, \quad A_z = A_z'=0.

From all these relations, the electromagnetic field can be computed as usual, through the definitions of potentials :math:`\mathbf{E}=-\nabla\Phi-\partial_t\mathbf{A}`, :math:`\mathbf{B}=-\nabla\times\mathbf{A}`:

.. math::
  \begin{eqnarray}
  E_x &=& -\partial_x \Phi - \partial_t A_x = -\partial_x \Phi + \beta_0^2 \partial_x \Phi = -\frac{1}{\gamma_0^2}\partial_x \Phi,\\ 
  E_y &=& -\partial_y \Phi - \partial_t A_y = -\partial_y \Phi,\\ 
  E_z &=& -\partial_z \Phi - \partial_t A_z = -\partial_z \Phi,\newline\\
  B_x &=& \partial_y A_z - \partial_z A_y = 0 ,\\ 
  B_y &=& \partial_z A_x - \partial_x A_z = \partial_z A_x = \beta_0 \partial_z \Phi = - \beta_0 E_z,\\   
  B_z &=& \partial_x A_y - \partial_y A_x = - \partial_y A_x = - \beta_0 \partial_y \Phi = \beta_0 E_y,
  \end{eqnarray} 

or in more compact form: :math:`\mathbf{E}=\left( -\frac{1}{\gamma_0^2}\partial_x \Phi, -\partial_y \Phi,-\partial_z \Phi \right)`, :math:`\mathbf{B}=\beta_0\mathbf{\hat{x}}\times\mathbf{E}`. 
  
From the previous equations, it can be inferred that, in a 1D cartesian geometry, the fields computed through this procedure equal those obtained through the standard Poisson's problem. 
This can also be inferred from the relativistic transformations of fields, which conserve the :math:`x` components of the electromagnetic fields for boosts in the :math:`x` direction. 

In the case of azimuthal Fourier decomposition (``geometry = "AMcylindrical"`` in the namelist - see :doc:`namelist`, :doc:`algorithms`), given the linearity of the relativistic Poisson equation, the full equation
can be decomposed in azimuthal modes, with the correspondent mode component of the charge density :math:`-\tilde{\rho}^m` as source term.

The relativistic Poisson equation for the potential component :math:`\tilde{\Phi}^m` of the mode :math:`m` in this  geometry is thus:

.. math::
  :label: RelPoissonModes

  \left[ \frac{1}{\gamma^2_0}\partial^2_x\tilde{\Phi}^m+\frac{1}{r}\partial_r\left(r\partial_r\tilde{\Phi}^m\right)-\frac{m^2}{r^2}\tilde{\Phi}^m \right]\Phi = -\tilde{\rho}^m.

Solving each of these relativistic Poisson's equations allows to initialize the azimuthal components of the electromagnetic fields:

.. math::
  \begin{eqnarray}
  \tilde{E}^m_x &=& -\frac{1}{\gamma_0^2}\partial_x \tilde{\Phi}^m,\\ 
  \tilde{E}^m_r &=& -\partial_r \tilde{\Phi}^m, \\ 
  \tilde{E}^m_{\theta} &=& \frac{im}{r} \tilde{\Phi}^m,\newline\\
  \tilde{\mathbf{B}}^m &=& \beta_0\mathbf{\hat{x}}\times\tilde{\mathbf{E}}^m.
  \end{eqnarray} 


----

Recommendations for relativistic species field initialization
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

In :program:`Smilei`, each species can independently benefit from this field initialization procedure. Its field will be initialized when the species will start to move, in order not to interfere with the other species' dynamics. 
The initialized fields will be superimposed to the electromagnetic fields already present in the simulation. To have physically meaningful results, we recommend to place a species which requires this method of field initialization far from other species, otherwise the latter could experience instantaneous unphysical forces by the relativistic species’ fields.

Remember that the transverse field of a moving charge with relativistic factor :math:`\gamma` is greater than the electrostatic transverse field of that charge, by a factor :math:`\gamma`. 
This means that for highly relativistic particles, you will need to use a transversely large simulation window to let the field decrease enough to reduce border effects during its propagation. 
A complete absence of boundary effects in this case would be provided by perfectly absorbing boundary conditions, which are not implemented yet in the code. 
If the relativistic species propagates in a plasma, these border effects could be partially screened by the plasma.

A relativistic mean velocity in the :math:`x` direction and a negligible energy spread are assumed in the hypotheses of this procedure, so the user must ensure these conditions when defining the species requiring field initialization in the namelist. 
The procedure could be extended to non-monoenergetic species, dividing the species particles in monoenergetic energy bins and then superimposing the fields by each of the monoenergetic bins, computed with the same procedure. 
At the moment, this energy binning technique is not available in :program:`Smilei`.  



----

References
^^^^^^^^^^

.. [Vay2008] `J.-L. Vay, Physics of Plasmas 15, 056701 (2008) <https://doi.org/10.1063/1.2837054>`_

.. [Londrillo2014] `P. Londrillo, C. Gatti and M. Ferrario, Nucl. Instr. and Meth. A 740, 236-241 (2014) <https://doi.org/10.1016/j.nima.2013.10.028>`_

.. [Massimo2016] `F. Massimo, A. Marocchino and A. R. Rossi, Nucl. Instr. and Meth. A 829, 378-382 (2016) <https://doi.org/10.1016/j.nima.2016.02.043>`_

.. [Marocchino2018] `A. Marocchino, E. Chiadroni, M. Ferrario, F. Mira and A.R. Rossi, Nucl. Instr. and Meth. A (2018) <https://doi.org/10.1016/j.nima.2018.02.068>`_




