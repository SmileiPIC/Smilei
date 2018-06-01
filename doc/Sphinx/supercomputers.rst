Referenced supercomputers
-------------------------

As it has been mentionned in :doc:`Install <installation>`, some compilation informations
(compiler options, environment) can be stored in dedicated files, stored in ``scripts/CompileTools/machine``.

It is especially the case for some supercomputers :

.. code-block:: bash
                
   [~/smilei]$ ls scripts/CompileTools/machine
   curie
   irene
   jureca
   occigen
   ...

Those files are composed of commands which must be executed by the user, then by compiler options automatically added to the default options of the makefile.

Commands are printed using the ``make machine=`` **target** ``help`` command :

.. code-block:: bash
                   
   [~/smilei]$ make machine=occigen help
   ...
   Machine comments for occigen:
   # module purge
   # module load intel intelmpi hdf5/1.8.18 qt/4.8.6 python/2.7.12 mesa/17.2.4 VTK/7.0.0

                
