.. _tablePage:

Generation of the external tables
--------------------------------------------------------------------------------

Several physical mechanisms can use external tables to work:

* Radiation loss and photon emission via nonlinear inverse Compton scattering (see :doc:`radiation_loss`)
* Electron-positon pair creation via the Multiphoton Breit-Wheeler (see :doc:`multiphoton_Breit_Wheeler`)

An external tool called :program:`smilei_tables` is available to generate these tables.

Installation
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

The C++ sources of this tool is located in `tools/tables`.

Required dependencies are the following:

* A C++11 compiler
* A MPI library
* HDF5 installed at least in serial
* Boost

Boost is a C++ library that provides efficient advanced mathematical functions.
This is the only dependency not required to install :program:`Smilei`.
This library can be easily installed manually on Linux, MacOS or Windows systems.
It is also available via different package managers (Debian, Homebrew).
The environment variable `BOOST_ROOT` must be defined.

The tool can be then installed using the makefile and the argument `tables`:

.. code-block:: bash

  make tables

The compilation generates an executable called :program:`smilei_tables` on the root of the repository.

Execution
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

The tool works with command line arguments.
For each physical mechanism, :program:`smilei_tables` generates all tables for this mechanism.
The first argument therefore corresonds to the physical mechanism:

* Nonlinear inverse Compton scattering: `nics`
* Multiphoton Breit-Wheeler: `mbw`
* For help: `-h` or `--help`

.. code-block:: bash

  mpirun -np <number of processes> ./smilei_tables -h

Then, once the physical mechanism selected, the following arguments are numerical parameters for the table generation.
For each physical argument, `-h` or `--help` gives the full list of arguments.

**For Nonlinear inverse Compton Scattering:**

.. code-block:: bash

  mpirun -np <number of processes> ./smilei_tables nics -h

  _______________________________________________________________________

  Smilei Tables
  _______________________________________________________________________

  You have selected the creation of tables for the nonlinear inverse Compton scattering.

  Help page specific to the nonlinear inverse Compton Scattering:

  List of available commands:
  -h, --help                       print a help message and exit.
  -s, --size       int int         respective size of the particle and photon chi axis. (default 128 128)
  -b, --boundaries double double   min and max of the particle chi axis. (default 1e-3 1e3)
  -e, --error      int             compute error due to discretization and use the provided int as a number of draws. (default 0)
  -t, --threshold  double          Minimum targeted value of xi in the computation the minimum particle quantum parameter. (default 1e-3)
  -p, --power      int             Maximum decrease in order of magnitude for the search for the minimum particle quantum parameter. (default 4)
  -v, --verbose                    Dump the tables

**For multiphoton Breit-Wheeler:**

.. code-block:: bash

  mpirun -np <number of processes> ./smilei_tables mbw -h

  _______________________________________________________________________

  Smilei Tables
  _______________________________________________________________________

  You have selected the creation of tables for the multiphoton Breit Wheeler process.

  Help page specific to the multiphoton Breit-Wheeler:

  List of available commands:
  -h, --help                       print a help message and exit.
  -s, --size       int int         respective size of the photon and particle chi axis. (default 128 128)
  -b, --boundaries double double   min and max of the photon chi axis. (default 1e-2 1e2)
  -e, --error      int             compute error due to discretization and use the provided int as a number of draws. (default 0)
  -t, --threshold  double          Minimum targeted value of xi in the computation the minimum photon quantum parameter. (default 1e-3)
  -p, --power      int             Maximum decrease in order of magnitude for the search for the minimum photon quantum parameter. (default 4)
  -v, --verbose                    Dump the tables

The tables are generated where the code is executed using HDF5 with the following names:

* Nonlinear inverse Compton Scattering: `radiation_tables.h5`
* multiphoton Breit-Wheeler: `multiphoton_breit_wheeler_tables.h5`

Precomputed tables
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

We have computed some tables with several levels of discretizations that you can download here.

256 points
"""""""""""

.. code-block:: bash

  mpirun -np <number of processes> ./smilei_tables nics -s 256 256 -b 1e-4 1e3
  
radiation_tables.h5

.. code-block:: bash

  mpirun -np <number of processes> ./smilei_tables mbw -s 256 256 -b 1e-2 1e2

multiphoton_breit_wheeler.h5

512 points
"""""""""""

.. code-block:: bash

  mpirun -np <number of processes> ./smilei_tables nics -s 256 256 -b 1e-4 1e3
  
`radiation_tables.h5 <http://mdls-internet.extra.cea.fr/projects/Smilei/uploads/tables_512/radiation_tables.h5>`_

.. code-block:: bash

  mpirun -np <number of processes> ./smilei_tables mbw -s 512 512 -b 1e-2 1e2

`multiphoton_breit_wheeler_tables.h5 <http://mdls-internet.extra.cea.fr/projects/Smilei/uploads/tables_512/multiphoton_breit_wheeler_tables.h5>`_

Python visualization scripts
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Detailed description of the tables
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Nonlinear Inverse Compton Scattering
""""""""""""""""""""""""""""""""""""

The file `radiation_tables.h5` is used for the nonlinear inverse Compton scattering radiation
mechanism described in :doc:`the dedicated section <radiation_loss>`.
It contains the `integfochi` table that represents
the integration of the synchortron emissivity of Ritus *et al*:

.. math::
  :label: eq_integfochi

  \int_{0}^{\chi_\pm} \frac{2 x}{ 3 \chi_\pm^2} \left[ \int_{2y}^{+\infty}{K_{1/3(y)}dy} - \frac{2 + 3 x y}{2} K_{2/3}(\nu) \right] dx

where

.. math::
  :label: eq_y
  
  y = \frac{x}{3 \chi_\pm (\chi_\pm - x)}

It is used by the Monte-Carlo method to compute the radiation emission cross-section.

.. _nics_integration_F_over_chi:

.. figure:: _static/nics/nics_integration_F_over_chi.png
  :width: 15cm

  Plot of the integfochi table for a particle quantum parameter ranging
  from :math:`\chi = 10^{-4}` to :math:`10^{3}` using the pre-computed table of 512 points.
