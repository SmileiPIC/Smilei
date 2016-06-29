Releases
--------

The code is hosted in two location:

1. All major releases are now available as compressed archive on `Github <https://github.com/SmileiPIC/Smilei/releases>`_.
  
  You're welcome to fork, develop and pull request as described `here <https://help.github.com/articles/fork-a-repo/>`_.

2. The developper *bleeding edge* official *Git* repository is `Gitlab <https://llrgit.in2p3.fr/smilei/smilei>`_.

  This repository is private and you will need to ask for a password first, but the process is fast and straightforward.


----

Release 2.2
^^^^^^^^^^^

This release includes:

* Patched-based parallelisation strategy and dynamic load balancing.
* Diagnostics have been strongly improved to deal at best with this parallelisation strategy.
* Improved readability with the new python input file format: all SMILEI variables are encapsulatedin a `Main` block (you will thus have to change your input file accordingly).
* SMILEI is now also shared on GitHub. Gitlab at LLR remains as main development platform.


----

Release 2.0
^^^^^^^^^^^

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


