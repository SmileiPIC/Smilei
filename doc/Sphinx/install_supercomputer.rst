Installation on supercomputers
-------------------------------

On a large cluster, refer to the administrator to install the requirements
and to choose the compilation options.

For a few machines, we have compiled the instructions in dedicated
files in the folder ``scripts/CompileTools/machine``:

- Cori, Haswell partition: `cori_hsw`
- Cori, Knights Landing partition: `cori_knl`
- Curie: `curie`
- Frioul: `frioul`
- Irene Joliot-Curie, KNL partition: `joliot_curie_knl`
- Irene Joliot-Curie, Skylake partition: `joliot_curie_skl`
- Jean Zay, Cascadelake: `jean_zay`
- Jureca, Haswell: `jureca`
- Marconi, Broadwell: `marconi_bdw`
- Marconi, Knights Landing: `marconi_knl`
- Occigen, Haswell: `occigen`
- Poincare, Sandy Bridge: `poincare`
- Tornado, Broadwell: `tornado`

These flags contain compiler flags specific for these machines to make Smilei the most efficient as possible.
We have also compiled instructions for some common architectures that you can use if your machine uses one of them:

- Intel Cascadelake processors: `cascadelake`
- Intel Skylake processors: `skylake`
- Intel Knights Landing processors: `knl`
- Intel Broadwell processors: `broadwell`
- Intel Broadwell processors: `haswell`

These files are composed of commands, which must be executed by the user,
followed by compiler options automatically added to the default options of the makefile.

To print out the commands to be executed, type ``make machine=target help``.
See, for instance:

.. code-block:: bash
                   
   $ make machine=occigen help
   ...
   Machine comments for occigen:
   # module purge
   # module load intel intelmpi hdf5/1.8.18 qt/4.8.6 python/2.7.12 mesa/17.2.4 VTK/7.0.0

After copying and pasting those commands to the terminal, you can use the command
``make machine=target`` to compile Smilei. For instance:

.. code-block:: bash

  $ make machine=occigen


If your machine is not referenced in this list, you must contact your administrator
for help on the installation. You may submit your installation instructions
to the Smilei repository in order to add your machine to the list.
