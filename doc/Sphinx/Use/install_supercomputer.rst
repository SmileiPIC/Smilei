Installation on supercomputers
-------------------------------

On a large cluster, refer to the administrator to install the requirements
and to choose the compilation options.

For a few existing machines, we provide instructions in the folder
``scripts/compile_tools/machine``. Each file contains compiler flags
and environment setup to optimize Smilei's performance.

.. rst-class:: fancy

+---------------------------------------------------------------------------+-----------------------------------------------------------------------+
| `Archer2 <https://www.archer2.ac.uk/>`_                                   | | Archer2 (GNU compiler): ``archer2``                                 |
+---------------------------------------------------------------------------+-----------------------------------------------------------------------+
| `Cori <http://docs.nersc.gov/systems/cori>`_                              | | Haswell: ``cori_hsw``                                               |
|                                                                           | | KNL: ``cori_knl``                                                   |
+---------------------------------------------------------------------------+-----------------------------------------------------------------------+
| `Frioul <http://frioul.int.univ-amu.fr>`_                                 | | ``frioul``                                                          |
+---------------------------------------------------------------------------+-----------------------------------------------------------------------+
| `FUGAKU <https://www.fugaku.r-ccs.riken.jp>`_                             | | Fujitsu compiler in trad mode :  ``fugaku_fujitsu_tm``              |
|                                                                           | | Fujitsu compiler in clang mode :  ``fugaku_fujitsu_cm``             |
+---------------------------------------------------------------------------+-----------------------------------------------------------------------+
| `Joliot-Curie <http://www-hpc.cea.fr/en/complexe/tgcc-JoliotCurie.htm>`_  | | KNL (Intel compiler): ``joliot_curie_knl``                          |
|                                                                           | | Skylake (Intel compiler): ``joliot_curie_skl``                      |
|                                                                           | | Rome (Intel compiler): ``joliot_curie_rome``                        |
|                                                                           | | A64FX with the GNU compiler: ``joliot_curie_gnu_a64fx``             |
|                                                                           | | A64FX with the ARM compiler: ``joliot_curie_arm_a64fx``             |
|                                                                           | | A64FX with the Fujitsu compiler: ``joliot_curie_fujitsu_a64fx``     |
+---------------------------------------------------------------------------+-----------------------------------------------------------------------+
| `Jean Zay <http://www.idris.fr/jean-zay>`_                                | | Cascadelake: ``jean_zay``                                           |
+---------------------------------------------------------------------------+-----------------------------------------------------------------------+
| `Jureca <http://apps.fz-juelich.de/jsc/hps/jureca>`_                      | | Haswell: ``jureca``                                                 |
+---------------------------------------------------------------------------+-----------------------------------------------------------------------+
| `Marconi <http://www.hpc.cineca.it/hardware/marconi>`_                    | | Broadwell: ``marconi_bdw``                                          |
|                                                                           | | KNL: ``marconi_knl``                                                |
+---------------------------------------------------------------------------+-----------------------------------------------------------------------+
| `Occigen <http://www.cines.fr/calcul/materiels/occigen>`_                 | | Haswell: ``occigen``                                                |
+---------------------------------------------------------------------------+-----------------------------------------------------------------------+
| `Ruche <https://mesocentre.pages.centralesupelec.fr/user_doc/>`_          | | Cascadelake (Intel): ``ruche``                                      |
+---------------------------------------------------------------------------+-----------------------------------------------------------------------+
| `Stampede <https://www.tacc.utexas.edu/systems/stampede>`_                | | KNL: ``stampede2_knl``                                              |
|                                                                           | | skylake: ``stampede2_skylake``                                      |
+---------------------------------------------------------------------------+-----------------------------------------------------------------------+


We also provide instructions for some common architectures:

- Intel Cascadelake processors: ``cascadelake``
- Intel Skylake processors: ``skylake``
- Intel Knights Landing processors: ``knl``
- Intel Broadwell processors: ``broadwell``
- Intel Haswell processors: ``haswell``

All these files contain:

* Commented commands that must be executed manually by the user
* Compiler options automatically accounted for during compilation

To print out the commands to be executed, type ``make machine=target help``.
See, for instance:

.. code-block:: bash

   $ make machine=occigen help
   ...
   Machine comments for occigen:
   # module purge
   # module load intel intelmpi hdf5/1.8.18 qt/4.8.6 python/2.7.12 mesa/17.2.4 VTK/7.0.0

After copying and pasting those commands to the terminal, you can use the
command ``make machine=target`` to compile Smilei. For instance:

.. code-block:: bash

  $ make machine=occigen


If your machine is not in this list, please contact your administrator
for help on the installation. You may submit your installation instructions
to the Smilei repository so that we can add your machine to the list.
