Troubleshoot
---------------------

If you encounter issues running your simulation, you can ask for help in the chat room on Element or publish an ``Issue`` on `GitHub <https://github.com/SmileiPIC/Smilei/issues>`_.
From previous users experience, the cause of most of the issues can be found performing some basic checks, listed in the following. 
Also, performing these simple tests helps us while finding an answer to your ``Issue``.

----

Simulation not starting / Error running a simulation
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Your simulation does not start or has an error while running.

* Check that the ``smilei`` executable is called correctly, e.g. it is present in the input namelist folder or correctly called from there. See also :doc:`run`.
* Check that ``smilei`` runs correctly with at least one input namelist in the ``benchmarks`` folder.
* Use the ``smilei_test`` executable on your input namelist. Does it display error or warning messages?
* Check your :doc:`profiles` in the input namelist, e.g. Avoid ``NaN`` values if you read them from a file or generate them through a Python function.
* Try ``git pull``, ``make clean`` and compile again. Also, be sure that uou have already used ``make happi`` for postprocessing. All these operations must be performed in the installation folder on your machine to run the simulations.
* Change the numerical parameters in the input namelist e.g. resolution, timestep, particles per cell, number of patches. Does the error occur again? See the dedicated sections in :doc:`namelist`.
* Change the number of MPI process, OpenMP threads. Do the simulations in the ``benchmarks`` folder run correctly with only one MPI process/OpenMP thread?
* Check that your MPI and OpenMP configuration is correct through parallel Hello World tests or other basic parallel programs.
* Try running a reduced simulation (less particles per cell, coarser resolution).
* Try using the same input namelist, and/or a reduced version (less particles per cell, coarser resolution) on a different machine.
* If the simulation stops after its start, does the error occur always at the same iteration? With ``print_every`` in the input namelist you can change the iteration print frequency in the log file. 

----

New simulation does not run
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

You have already run successfully other similar simulations, but a new one gives an error at the start.

* Use the input namelist of the simulation that works correctly and then progressively change the input to arrive to the simulation you want to run. At each step check if something goes wrong and what is the change in the input namelist that caused the issue.
* Try running the new simulation with the ``smilei`` executable used for the simulation that worked correctly, if different. Changing the code version creates the problem?

----

Postprocessing error
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

You can run the simulation but you cannot open/read the results.

* Do ``git pull`` and ``make happi`` in your installation folder to have the postprocessing library ``happi`` ready for use. Afterwards, remember to close and reopen the Python interface you are using, e.g. IPython. See also :doc:`installation`.
* Check if you can correctly open the results of a simulation using one input namelist in the ``benchmarks`` folder.
* Carefully read the doc about the :doc:`post-processing` method you are trying to use.


----

Physical error in the results.
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

The physical results are not the ones you expect. 

* Read the doc on the physical methods you are using (e.g. :doc:`/Understand/collisions`,
  :doc:`/Understand/ionization`, :doc:`/Understand/laser_envelope`, ...).
  Are the underlying physical assumptions satisfied?
* Check that the units given in the input namelist are properly normalized.
  See also :doc:`/Understand/units`.
* Some physical processes like :doc:`/Understand/collisions`, :doc:`/Understand/ionization`
  need a reference frequency in SI in the ``Main`` block of the input namelist. Did you provide it?
  See also :doc:`namelist`.
* Check the CFL condition in the input namelist. See :doc:`/Understand/algorithms`
* See with the Scalar diagnostics (See :doc:`post-processing` ) if the kinetic energy ``Ukin``
  or electromagnetic energy ``Uelm`` display strange behaviours (e.g. exponential growths).
* Verify the overall consistency of the physical set-up, e.g. only immobile or almost immobile
  particles while using a Poisson solver.
* Verify that the physical initialization is correct. Should you use a classical or
  relativistic Poisson solver (See :doc:`/Understand/relativistic_fields_initialization`)
  for the initial fields. Is it necessary to use a Poisson solver?
* Check the presence of numerical effects running the simulation with different numerical
  parameters, e.g. changing the resolution, timestep, in the input namelist.
* If using the ``AMcylindrical`` geometry, check that the origin of the axes you are using
  in the input namelist is the one described in See :doc:`/Understand/azimuthal_modes_decomposition`.

----

Performances issues
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

The simulation is very slow / the performances are not as expected.

* Change the number of MPI process and OpenMP threads.
* Change the number of patches and/or their distribution in each direction.
  See also :doc:`/Understand/parallelization`.
* Check that ``LoadBalancing`` is activated in the :doc:`namelist`
  (if the physical set-up is suitable for its use). See also :doc:`/Understand/parallelization`.
* If using :doc:`/Understand/vectorization`, check that the compilation flags for vectorization
  were correctly used. See also :doc:`installation`.


