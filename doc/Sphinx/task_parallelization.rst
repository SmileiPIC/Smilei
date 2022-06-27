Task Parallelization
----------------------

Task parallelization is a method to spread the computing workload on the many cores
of a computer. Instead of splitting the *data* accross cores and apply the same task
to all these pieces, different *tasks* are split accross the cores sharing the same
data. This approach can often make the computation faster, especially
with non-uniform plasma distributions.

Task parallelization of macro-particle operations in Smilei (using OpenMP) is
published in [Massimo2022]_.

----

Motivation
^^^^^^^^^^

Usually, most of the computing time is spent on macro-particle operations.
Consequently, a non-uniform plasma distribution results in load imbalance:
some cpu cores are loaded with less macro-particles, and idly wait for the
other cores to finish their work.

Worse: adding more cpu cores will not result in significant speedup, as only
a few of them are performing most of the work. This "strong scaling" curve
(speed-up vs number of cpu-cores) starts to saturate. Several methods for
parallelism can increase the number of computing units where
the saturation occurs.

In :program:`Smilei`, by default, the data is split in *patches* (see
:doc:`parallelization`) and, when the environment variable ``OMP_SCHEDULE``
is set to ``dynamic``, the OpenMP scheduler dynamically assigns each patch
to each core. This provides for some load balancing at the MPI level, as
cores can work asynchronously on different patches.

This strategy implies that only 1 OpenMP thread can work on a given patch,
which includes potentially several ``Species`` and all the PIC operators
(interpolation, push, etc). These constraints can considerably slow down the
simulation in some situations (many species with non-uniform distribution,
and/or low number of patches per core).

A first solution is to split the data to a finer level: separate the
treatment of species, and split the patch in smaller structures (in Smilei,
patches are divided in ``clusters`` along the dimension ``x``). This
can improve the strong scaling results, but some constructs cannot be
parallelized with this data splitting (e.g. irregularly nested loops, recursion,
etc). The task parallelism has been introduced to answer these issues.

----

Task approach
^^^^^^^^^^^^^^^^^^^^^^^^^^^^

:program:`Smilei` exploits the task parallelization available since OpenMP 4.5.
The main idea is to split the work in smaller units that can be run asynchronously.

In addition to separated species treatment and patches split in clusters
(see :py:data:`cluster_width`), the macro-particle operators (interpolation, push, etc)
are defined as tasks. All the combinations of [operator-cluster-species-patch]
correspond to different tasks that can be run in parallel. 

As some tasks depend on other tasks, the dependency tree is provided to OpenMP so
that the tasks are dynamically assigned to OpenMP threads, in the correct order
(preventing race conditions). This is described in [Massimo2022]_.

----

Performance Results
^^^^^^^^^^^^^^^^^^^^^

Some results from [Massimo2022]_ are shown in the following.

A 2D uniform thermal plasma case show that with uniform macro-particle 
distributions the task-parallelization in :program:`Smilei` does not have a 
performance advantage.

.. _uniform_plasma:

.. figure:: _static/Scan_Uniform_Plasma_2D.png
    :width: 40%
    :align: center

    Performances with and without task parallelization in a uniform plasma case.

However, a 2D radiation pressure acceleration is an example of non-uniform 
macro-particle distribution where the task parallelization yields and advantage.

.. _radiation_pressure_rho:

.. figure:: _static/Radiation_Pressure_Rho.png
    :width: 50%
    :align: center

    Electron density divided by the critical density in a 2D radiation pressure 
    benchmark at 0 (left) and 1500 iterations (right). The non-uniformity of the 
    macro-particle distribution is present since the start of the simulation.

.. _radiation_pressure_perf:

.. figure:: _static/Scan_Radiation_Pressure_2D.png
    :width: 40%
    :align: center

    Performances with and without task parallelization in a 2D radiation 
    pressure acceleration case.

The scheduling of macro-particle operations without and with task parallelization
can be seen in the following figures.
Note how in the first Figure (without task parallelization), the end of the 
treatment of macro-particle operators (around 0.1 s) is determined by the 
OpenMP thread 0 of the MPI process 0. In the second Figure (with task parallelization),
the OpemMP thread 2 of MPI process 0 determines the end of the 
treatment of macro-particle operators (around 0.07 s). In this case, the finer 
decomposition given by the bins and the relaxation of the constraints involved
in the assignment of macro-particle operations to threads yields a shorter time
to the result.

.. _task_tracing_tasks_off:

.. figure:: _static/Radiation_pressure_develop_tracing.png
    :width: 50%
    :align: center

    Scheduling of macro-particle operations for the 2D radiation pressure benchmark, 
    4 MPI processes and 4 OpenMP threads, during iteration 1200,
    without task parallelization.

.. _task_tracing_tasks_on:

.. figure:: _static/Radiation_pressure_task_tracing.png
    :width: 50%
    :align: center

    Scheduling of macro-particle operations for the 2D radiation pressure benchmark, 
    4 MPI processes and 4 OpenMP threads, 4 bins per patch, during iteration 1200, 
    with task parallelization. The horizontal axis has been extended to the same 
    maximum value of the horizontal axis of the previous Figure to facilitate 
    the comparison.
