Perfectly Matched Layers
-------------------------------------

Perfectly Matched Layers (PML) are open boundary conditions for electromagnetic and laser envelope fields.
This means that a wave propagating through the boundary will be absorbed and not reintroduced into the simulation.
In that regard, it plays the same role as the "Silver-Muller" boundary conditions which is supported for electromagnetic fields only (not envelope).

At the cost of having a slightly higher computing cost than Silver-Muller, PML have several very interesting advantages:

* All resolved frequencies are equally absorbed.
* Absorption is efficient for any incident angle.
* The overall absorption of exiting waves is, in most cases, much better than Silver-Muller.

This leads to a more accurate result of course but can also relax constraints on the size of the simulation domain.
Very large transverse sizes to let the laser diffract more before it hits the boundary is no longer necessary for instance.

.. rubric:: Basic use

A PML is in fact several layers of absorbing media which are contiguous to the simulation domain and through which exiting waves are going.
The waves are progressively absorbed as they travel across this absorbing media.
The thickness of this media can be tuned with the parameter `number_of_pml_cells` and is in fact equal to 

.. math::
  :label: PMLthickness

  PML_{\rm thickness} = {\rm number\_of\_pml\_cells[0]}\times {\rm dx}, 

in the `x` direction.

Keep in mind that **the absorbing power of a PML increases with its physical thickness**.

This means that the more layers are set, the more absorbing the PML is.
But adding layers also increases computing time of course because it is similar to increasing the simulation domain size.

It also means that the absorption of the PML depends on the resolution `dx`. 
Higher resolution simulations will need more cells in order to keep the same thickness and therefore the same absorption as a less resolved one.

The default number of pml cells is 6.

.. rubric:: Rules of thumb for a good use of the PML

First of all, it is important to realize that open boundary conditions are only good at letting waves out of the domain.
It does not relax any constraint on the quality of the initial conditions and will not prevent reflections if these are physical.

The typical case is a vacuum/plasma interface at the boundary of the simulation box when a `remove` boundary condition is chosen for particles.
In that configuration, if the plasma density is high enough at the interface, there can be perfectly physical reflections that PML are not designed to mitigate.
We therefore advise to use plasma profiles whith decreasing densities close to the boundaries.

.. rubric:: Advanced settings for Cartesian Geometry

In the PML medium, the complex susceptibility is given by

.. math::
  :label: PMLsusceptibility

  s = \kappa + \frac{\sigma}{i\omega\epsilon_0}, 

where :math:`\kappa` and :math:`\sigma` can be chosen by the user.
They are functions of a single space variable which describes the normal position in the PML.
This variable ranges from 0 (the inner bound of the PML) to 1, the outer bound of the PML.

One profile per simulation dimension is required and the same profile is applied to both sides of a given dimension.

.. rubric:: Expert settings for AM geometry

In the specific case of the AM geometry, using arbitrary susceptibility profiles is trickier because an integration along the radial direction is necessary.
Therefore it is important that :program:`Smilei` knows not only the profiles of `sigma` and `kappa` but also a primitive of these profiles along the radial dimension.

It is up to the user to provide such primitives in the namelist using the following syntax (here is an example for `sigma`):

  | **Syntax 1:** ``[sigma, integrate_sigma]``, identical for all dimensions.
  | **Syntax 2:** ``[sigma_x, sigma_r, integrate_sigma_r]``,  different on each dimension.

The default profiles identical for all dimensions and are given by:

.. code-block:: python

    def sigma(u):
        return 20. * u**2  
    def integrate_sigma(u):
        return 20./3. * u**3  
    def kappa(u):
        return 1 + 79. * u**4  
    def integrate_kappa(u):
        return u + 79./5. * u**5  

.. rubric:: PML for the envelope model

For stability purposes, the PML boundaries for the envelope use frequency shifting which prevents from using arbitrary profiles for the susceptibility.
Therefore it is not possible to tune its profile.

Details of the method will be published soon.
