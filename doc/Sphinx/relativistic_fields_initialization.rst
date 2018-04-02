.. _relativisticfieldsinitializationPage:

Relativistic electromagnetic field initialization
--------------------------------------------------------------------------------
From Maxwell's equations it can be shown that, provided that continuity equation holds, i.e.

.. math::

  \partial_t \rho + \nabla \cdot \mathbf{J} = 0,

the quantities :math:`\nabla\cdot\mathbf{B}`, :math:`\nabla\cdot\mathbf{E}-\rho` do not change over time:

.. math::

  \partial_t \left( \nabla\cdot\mathbf{B} \right ) = 0, 

.. math::

  \partial_t \left( \nabla\cdot\mathbf{E}-\rho \right ) = \nabla\cdot\partial_t\mathbf{E}-\partial_t\rho = \nabla\cdot\left(\nabla\times\mathbf{B}-\mathbf{J}\right)-\partial_t\rho = - \left(\nabla\cdot\mathbf{J}+\partial_t \rho\right).

Thus, if a simulation starts with :math:`\rho\neq0`, the electromagnetic fields must be properly initialized. 

In the case of an static charge distribution, i.e. :math:`\rho\neq0`, :math:`\mathbf{J}=0`, the initial electrostatic potential :math:`\Phi` can be computed solving Poisson's equation:

.. math::

  \nabla^2 \Phi = -\rho,

and then be derived to find the initial electric field: :math:`\mathbf{E}=-\nabla\Phi`. The initial magnetic field :math:`\mathbf{B}` will be zero.

In general when the initial current :math:`\mathbf{J}` is not zero, the electric and magnetic fields will not be zero and the full set of Maxwell's equations or equivalently the potentials equations should be solved. 

For physical setups where a species is already relativistic when the simulation starts, e.g. a relativistic electron bunch, the initial electromagnetic fields can be computed through a simplified procedure, described in [Vay2008]_, [Londrillo2014]_ and [Massimo2016]_. 

An important assumption of this calculation is that the species is highly relativistic, moving in the positive :math:`x` direction, with negligible velocity spread. Under this hypothesis, the transverse components of the species current density are neglected and the four-current quadrivector can be written as:

.. math::

  \left(\mathbf{J},\rho c\right) = \left(\rho \beta_0 c, 0, 0, \rho c\right),

where :math:`\beta_0 c` is the initial mean velocity of the relativistic species. At least locally, the potentials :math:`\mathbf{A}`, :math:`\Phi` in the laboratory frame will be only function of :math:`x-\beta_0 c t`, as they are propagating with the species at uniform velocity.

In the relativistic species rest frame :math:`S'`, the charge distribution is static and the electrostatic potential in that reference frame :math:`\Phi'` is related to the charge density in that reference frame :math:`\rho'` through Poisson's equation:

.. math::
  :label: Poisson

  \nabla'^2 \Phi' = -\rho',

where the Laplacian operator is computed in the reference frame :math:`S'`:

.. math::
  
  \nabla'^2=\partial^2_{x'}+\partial^2_{y'}+\partial^2_{z'}`.

The vector potential in the species rest frame can be set to zero: :math:`\mathbf{A'}=0`. Through the above mentioned assumptions, it is possible to rewrite Eq. :eq:`Poisson` only in terms of laboratory frame quantities. 

Lorentz transformation of the four-vector :math:`\left(\mathbf{J},\rho c\right)` yields :math:`\rho'=\rho/\gamma_0`, where :math:`\gamma_0=1/\sqrt{1-\beta^2_0}` is the average Lorentz factor of the relativistic species. 
Similarly, the potential :math:`\Phi'` can be rewritten in terms of the potential in the laboratory frame: :math:`\Phi'=\Phi/\gamma_0`. The Lorentz back-transformation of coordinates

.. math::
  
  x=\gamma_0(x'+\beta_0 ct'),\quad  ct = \gamma_0(ct'+\beta x'), \quad y=y', \quad z=z'

allows to transform the derivatives in Eq. :eq:`Poisson` as 

.. math::
  
  \partial_{x'}=\gamma_0\left(\partial_x+\frac{\beta_0}{c}\partial_t\right), \quad \partial_{y'}=\partial_y, \quad \partial_{z'}=\partial_z. 

The partial derivative along the :math:`x'` direction can be further simplified, through the hypothesis of temporary dependence of all quantities on :math:`x-\beta_0 c t`, implying :math:`\partial_t=-\beta c\partial_x`:

.. math::
  
  \partial_{x'}=\frac{1}{\gamma_0}\partial_x. 

Equation :eq:`Poisson` can thus be rewritten as 

.. math::
  :label: RelPoisson

  \left( \frac{1}{\gamma^2_0}\partial^2_x+\nabla_{\perp}^2\right) \Phi = -\rho,

here referred to as relativistic Poisson's equation.

Once the potential :math:`\Phi` is found, we can compute all the components of the electromagnetic field, using again the relations :math:`\partial_t=-\beta c\partial_x`, :math:`\Phi'=-\Phi/\gamma_0` and the Lorentz back-transformation of the vector potential :math:`\mathbf{A}`:

.. math::
  
  A_x = \gamma_0(A_x'+\beta_0 \Phi'/c)=\gamma_0\beta_0 \Phi'/c=\beta_0\Phi/c,\quad A_y = A_y'=0, \quad A_z = A_z'=0,

.. math::

  E_x = -\partial_x \Phi - \partial_t A_x = -\partial_x \Phi + \beta_0^2 \partial_x \Phi = -\frac{1}{\gamma_0^2}\partial_x \Phi,\\ 

.. math::

  E_y = -\partial_y \Phi - \partial_t A_y = -\partial_y \Phi,\\ 

.. math::

  E_z = -\partial_z \Phi - \partial_t A_z = -\partial_z \Phi,\\

.. math::

  B_x = \partial_y A_z - \partial_z A_y = 0 ,\\ 

.. math::

  B_y = \partial_z A_x - \partial_x A_z = \partial_z A_x = \frac{\beta_0}{c} \partial_z \Phi = - \frac{\beta_0}{c} E_z,\\ 

.. math::
  
  B_z = \partial_x A_y - \partial_y A_x = - \partial_y A_x = - \frac{\beta_0}{c} \partial_y \Phi = \frac{\beta_0}{c} E_y.
   
  










----

References
^^^^^^^^^^

.. [Vay2008] `J.-L. Vay, Physics of Plasmas 15, 056701 (2008) <https://doi.org/10.1063/1.2837054>`_

.. [Londrillo2014] `P. Londrillo, C. Gatti and M. Ferrario, Nucl. Instr. and Meth. A 740, 236-241 (2014) <https://doi.org/10.1016/j.nima.2013.10.028>`_

.. [Massimo2016] `F. Massimo, A. Marocchino and A. R. Rossi, Nucl. Instr. and Meth. A 829, 378-382 (2016) <https://doi.org/10.1016/j.nima.2016.02.043>`_


