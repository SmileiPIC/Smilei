ParticleBinning units
-----------------------------------

When using happi to post-process the output from ``ParticleBinning`` diagnostics,
the raw data is converted to more convenient units.

In a first step, happi divides the raw data in each bin by the size of the bin.
For instance, if the :py:data:`axes` of the diagnostic are ``x`` and ``ekin``, the
raw data is divided by a bin size that has units :math:`L_r K_r`.

A second, more subtle conversion is done by happi.
The raw quantity stored in the output file has the units of the :py:data:`deposited_quantity`.
Let us take the most common example: ``deposited_quantity`` is a sum of
:ref:`macro-particle weights<Weights>`. As those weights
are not in units of density, but of density multiplied by hypervolume,
the raw data is in units of :math:`N_r L_r^D`, where :math:`D` is the
dimension of the simulation.
The correction applied by *happi* is to divides the data by an
hypervolume. The choice of this hypervolume depends on the :py:data:`axes`
of the diagnostic: for each direction ``x``, ``y`` or ``z``, if this
direction was not already accounted for in the first step, then divide
by the length of the box in that direction.

To be clearer, let us take an example of a 3D simulation with a 
``ParticleBinning`` diagnostic that has axes ``x`` and ``ekin``,
and ``deposited_quantity="weight"``. Happi will first divide the raw
data by the bin size, which can be written :math:`\Delta x \Delta K`.
Then, it divides by the hypervolume, which is :math:`L_y L_z` (the
length of the box along :math:`x` is already accounted for when dividing
by the bin size). Finally, the units of the results will be
:math:`N_r L_r^3/(L_r K_r L_r^2) = N_r/K_r`: a density per unit energy.



