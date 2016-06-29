Releases
--------

All major releases are now available on `Github <https://github.com/SmileiPIC/Smilei>`_. 

If you are looking for an access to the *bleeding edge* official *Git* repository, check it out `here <https://llrgit.in2p3.fr/smilei/smilei>`_.
You will need to ask for a password first, but the process is fast and straightforward.


----

.. _latestVersion:

Current release 2.2
^^^^^^^^^^^^^^^^^^^

This release includes:

 * Patched-based parallelisation strategy and dynamic load balancing.
 * Diagnostics have been strongly improved to deal at best with this parallelisation strategy.
 * Improved readability with the new python input file format: all SMILEI variables are encapsulatedin a `Main` block (you will thus have to change your input file accordingly).
 * SMILEI is now also shared on GitHub. Gitlab at LLR remains as main development platform.


----

Release 2.0
^^^^^^^^^^^

**Download**: `Smilei v2.0 (pre-release) <_downloads/smilei-v2.0.tar.gz>`_

A second release is **in preparation** featuring great improvements:

* **state-of-the-art** dynamic load balancing
* full *python* namelist, allowing for complex, user-friendly input
* 3D cartesian geometry
* external fields and antennas
* binary collisions
* new diagnostics
* *python* scripts for post-processing
* various improvements on stability and error management

We greatly appreciate external users trying this code and giving feedback.
You can :doc:`contact us <partners>` to become a developer of the official releases.


----

Release 1.0
^^^^^^^^^^^

**Download**: `Smilei v1.0 <_downloads/smilei-v1.0.tar.gz>`_

After nearly two years of development, :program:`Smilei` v1.0 offers some nice features:

* 1D & 2D cartesian geometries
* moving-window
* hybrid MPI-OpenMP parallelization
* field ionization
* some python diagnostics
* open-source

.. warning::
  This version does not have associated documentation.
  Please refer to the examples.


