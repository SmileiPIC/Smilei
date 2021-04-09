Developer Zone
-----------------------------

Introduction
^^^^^^^^^^^^^^^^^^^^^^^^^^^

Smilei is a C++ code that uses relatively simple C++ features for modularity
and conveniency for non-advanced C++ users.

The repository is composed of the following directories:

- ``Licence``: contains code licence information
- ``doc``: conatins the Sphinx doc files
- ``src``: contains all source files
- ``happi``: contains the sources of the happi Python tool for visualization
- ``benchmarks``: contains the benchmarks used by the validation process. these becnhamrks are also examples for users.
- ``scripts``: contains multiple tool scripts for compilation and more
  - ``compile_tools``: contains scripts and machine files used by the makefile for compilation
- ``tools``: contains some additional programs for Smilei
- ``validation``: contains the python scripts used by the validation process

General concept and vocabulary
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

This section presents some implicit notions to understand the philosophy of the code.

Notion of data container
""""""""""""""""""""""""""""""

Data containers are classes (or sometime just structures) used to store a specific type of data, often considered as raw data such as particles or fields.
Some methods can be implemented in a data container for managing or accessing the data.

.. _dataContainer:

.. figure:: _static/figures/data_container.png
  :width: 5cm

  Data container.

Notion of operators
""""""""""""""""""""""""""""""

An operator is a class that operates on input data to provide a processed information.
Input data can be parameters and data containers.
Output data can be processed data from data containers or updated data containers.
An operator is a class functor (overloadind of the ``()`` ).
Sometime, operator provides additional methods called wrappers to provide differents simplified or adapted interfaces.
An operator do not store data or temporarely.
for instance, the particle interpolation, push and proection are operators.

.. _operator:

.. figure:: _static/figures/operator.png
  :width: 10cm

  Operator.

Notion of domain parts
""""""""""""""""""""""""""""""

Domain parts are classes that represents some specific levels of the domain decomposition.
They can be seen as high-level data container or container of data container.
They contain some methods to handle, manange and access the local data.
For instance, patches and ``species`` are domain parts:

- ``species`` contains the particles.
- ``patches`` contains ``species`` and fields.

Notion of factory
""""""""""""""""""""""""""""""

Some objects such as operators or data containers have sereral variations.
For this we use inheritance.
A base class is used for common parameters and methods and derived classes are used for all variations.
The factory uses user-defined input parameters to determine the right derive class to choose and initiate them.
For instance, there are several ``push`` operators implemented all derived from a base ``push`` class.
The ``push`` factory will determine the right one to use.

Other
""""""""""""""""""""""""""""""

Some classes are used for specific actions in the code such as the initilization process.

Domain decomposition and parallelism
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

The simulation domain is divided multiple times following a succession of decomposition levels.
The whole domain is the superimposition of different grids for each electromagnetic field components
and macro-particules.
Let us represent schematically the domain as an array of cells as in Fig. .
Each cell contains a certain population of particles.




Data structures
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Patches
""""""""""""""""""""""""""""""


The basic PIC loop implementation
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
