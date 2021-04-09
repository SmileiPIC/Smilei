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
  :width: 15cm

  Data container.

Notion of operators
""""""""""""""""""""""""""""""

An operator is a class that operates on input data to provide a processed information.
Input data can be parameters and data containers.
Output data can be processed data from data containers or updated data containers.
An operator is a class functor (overloadind of the `()` ).
Sometime, operator provides additional methods called wrappers to provide differents simplified or adapted interfaces.
An operator do not store data or temporarely.

Notion of domain parts
""""""""""""""""""""""""""""""

Domain parts are classes that represents some specific levels of the domain decomposition.
They can be seen as high-level data container or container of data container.
They contain some methods to handle, manange and access the local data.
for instance, patches and species are domain parts.

Notion of factories
""""""""""""""""""""""""""""""

Data structures
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^


The basic PIC loop implementation
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
