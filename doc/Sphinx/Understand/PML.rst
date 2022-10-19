PML: Perfectly Matched Layers
-------------------------------------

Open boundary conditions

.. rubric:: Parametrization in Cartesian Geometry

* number_of_pml_cells
* sigma
* kappa

.. rubric:: Rules of thumb for a good use of the PML

How to properly chose sigma and kappa ?

.. rubric:: PML in AM geometry

Need to provide a primitive for the radial profiles

.. rubric:: PML for the envelope model

For stability purposes, the PML boundaries for the envelope use frequency shifting which prevents from using arbitrary profiles. 
Therefore it is not possible to tune the profiles.

Details of the method are found here:
`in this paper <https://link_to_guillaume's paper>`_.
