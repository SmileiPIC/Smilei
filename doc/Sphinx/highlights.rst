Highlights
----------

Electron Acceleration
^^^^^^^^^^^^^^^^^^^^^

Below, an example of electron acceleration by laser wakefield.
The figure represents the evolution of the electronic density in time. 
A hotspot of electron is created behind the bubble.

.. raw:: html

    <video controls="controls">
    <source src="_static/Rho_electron1long.ogg" type="video/ogg" />
    </video>

----

Scalability in a wakefile acceleration simulation
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Wakefield-acceleration of electrons in an underdense plasma creates a
hotspot of electrons, which makes the simulation strongly imbalanced.
This spot represent a large part of the total calculations, so that
more computing power should be allocated on it.

Please refer to the doc :doc:`parallelization` to learn the basics of the
parallelization techniques employed in this section.

.. rubric :: 1. OpenMP

In a local area around this hotspot, OpenMP is able to manage the computing
resources to make the overall simulation faster. The following figure shows
the evolution of the time to calculate 100 iterations, as a function of time.
Each line corresponds to a different partition of the box in terms of
MPI processes and OpenMP threads: :math:`N\times M`, where :math:`N` is 
the total number of MPI processes, and :math:`M` is the number of threads
in each MPI process.

.. image:: _static/openMP_balancing.png
    :width: 500px

Using more OpenMP threads per MPI process (while keeping the total number
of threads constant) clearly reduces the simulation time, because the 
computing power is balanced within each MPI region.


.. rubric :: 2. Dynamic load balancing between MPI processes

At the global simulation scale, OpenMP cannot be used to smoothen the balance.
Instead, a dynamic load balancing (DLB) algorithm periodically exchanges pieces of 
the simulation box (*patches*) between MPI processes, so that each MPI
process owns a fair amount of the simulation load. The following figure
shows how this balancing reduces the time of the simulation.

.. image:: _static/DLB_balancing.png
    :width: 500px

The red curve is the best situation obtained in the previous section, while
the black curve corresponds to the DLB algorithm enabled.

The portion of the box belonging to each MPI process varies when the load balancing
occurs. The following figure shows how each of these portions evolve with time.

.. image:: _static/Patch_loadcomparision.jpg

The four panels correspond to four timesteps during the simulation.
The colorscale represents the log-scaled load of each patch.
The black lines show the borders of each MPI process' portion of the box.
The MPI processes that are close to the hotspot tend to handle a smaller portion
of the box.