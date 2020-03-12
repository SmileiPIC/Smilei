Installation on supercomputers
-------------------------------

On a large cluster, refer to the administrator to install the requirements
and to choose the compilation options.

For a few existing machines, we provide instructions in the folder
``scripts/CompileTools/machine``. Each file contains compiler flags
and environment setup to optimize Smilei's performance.

.. rst-class:: fancy

+--------------------------------------------------------------------+---------------------------------+
| `Cori <docs.nersc.gov/systems/cori>`_                              | | Haswell: ``cori_hsw``         |
|                                                                    | | KNL: ``cori_knl``             |
+--------------------------------------------------------------------+---------------------------------+ 
| `Frioul <frioul.int.univ-amu.fr>`_                                 | | ``frioul``                    |
+--------------------------------------------------------------------+---------------------------------+
| `Joliot-Curie <www-hpc.cea.fr/en/complexe/tgcc-JoliotCurie.htm>`_  | | KNL: ``joliot_curie_knl``     |
|                                                                    | | Skylake: ``joliot_curie_skl`` |
+--------------------------------------------------------------------+---------------------------------+
| `Jean Zay <www.idris.fr/jean-zay>`_                                | | Cascadelake: ``jean_zay``     |
+--------------------------------------------------------------------+---------------------------------+
| `Jureca <apps.fz-juelich.de/jsc/hps/jureca>`_                      | | Haswell: ``jureca``           |
+--------------------------------------------------------------------+---------------------------------+
| `Marconi <www.hpc.cineca.it/hardware/marconi>`_                    | | Broadwell: ``marconi_bdw``    |
|                                                                    | | KNL: ``marconi_knl``          |
+--------------------------------------------------------------------+---------------------------------+
| `Occigen <www.cines.fr/calcul/materiels/occigen>`_                 | | Haswell: ``occigen``          |
+--------------------------------------------------------------------+---------------------------------+


We also provide instructions for some common architectures:

- Intel Cascadelake processors: ``cascadelake``
- Intel Skylake processors: ``skylake``
- Intel Knights Landing processors: ``knl``
- Intel Broadwell processors: ``broadwell``
- Intel Broadwell processors: ``haswell``

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
