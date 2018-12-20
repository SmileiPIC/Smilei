Contribute
----------

Contributions to the development of :program:`Smilei` are welcome.

To report bugs, please create an issue on the
`GitHub website <https://github.com/SmileiPIC/Smilei/issues/new>`_.

To develop new features, you may freely *fork*
the source code from the same website. When
your modifications to the code are ready, you can make a `pull request
<https://github.com/SmileiPIC/Smilei/pulls>`_ for review and
merge with the main repository.

Guidelines for new developments are:

* Write clear, commented code, or self-explanatory.
* Write the documentation corresponding to the new features, if any.
* Make validation cases, and reference data, corresponding to the added features.

----

Write documentation
^^^^^^^^^^^^^^^^^^^

The documentation you are currently reading is written in the
`reStructuredText <www.sphinx-doc.org/en/stable/rest.html>`_ (rST) language, and included
in the main :program:`Smilei` repository. This is a fairly simple markup language. You
can see examples from the source files, which are located in the
``doc/Sphinx`` folder and have the extension ``.rst``.

To transform it into an *html* website, it is
processed using the `sphinx <www.sphinx-doc.org>`_ python package that you may have to
install.
If you have sphinx installed, you may simply go to the
main :program:`Smilei` folder from a command line terminal, then run the command

.. code-block:: bash

   make doc

This creates a local *html* website accessible in the ``build/html/`` folder. Simply
open the ``build/html/index.html`` file in your favorite web browser.

To document a new feature, please modify the file ``namelist.rst`` to indicate the
syntax changes in the input file. If the feature requires detailed physical or numerical
background, you may add a new page in the "Understand" section of the website.
To do that, create a new ``.rst`` file, then reference it in the table of contents
located in ``index.rst``.

----

.. _Validation:

Validate your changes
^^^^^^^^^^^^^^^^^^^^^

:program:`Smilei` has a system of test cases combined with reference results that must be
validated, ideally, for every *push* submitted to the git repository.
These **test cases** are located in the ``benchmarks/`` folder.

Each benchmark has an associated **validation file**, written in *python*, which contains
instructions on how to produce an analysis of the results that can validate that particular
benchmark. The validation files are located in the ``validation/analyses/`` folder.
They have the same name as the benchmarks, with the prefix ``validate_``.

Once a benchmark has been run, the corresponding ``validate_*`` file is run in *python*
to compare the analysis results with a **reference file** located in the folder
``validation/references/``. Note that the analysis produced by the ``validate_*`` file 
can also be used to generate the reference file the first time.

**When you code a new feature, you must provide a new benchmark, the corresponding
analysis and reference file**

To make this process easier, a *python* script is available.

.. rubric:: How do I use the ``validation.py`` script?

The script ``validation/validation.py`` can do three things:

* generate validation reference(s) for given benchmark(s)
* compare benchmark(s) to their reference(s)
* show visually differences between benchmark(s) and their reference(s)

Usage:

..
  
  .. code-block:: bash
  
    python validation.py [-c] [-h] [-v] [-o <nOMP>] [-m <nMPI>] [-b <bench> [-g | -s]] [-r <nRestarts>]
  
  * | Option ``-b <bench>``:  
    | ``<bench>`` : benchmark(s) to validate. Accepts wildcards.  
    | ``<bench>=?`` : ask input for a benchmark  
    | DEFAULT : All benchmarks are validated.   
  
  * | Option ``-o <nOMP>``:
    | ``<nOMP>`` : number of OpenMP threads used for the execution
    | DEFAULT : 4  
  
  * | Option ``-m <nMPI>``:
    | ``<nMPI>`` : number of MPI processes used for the execution
    | DEFAULT : 4
  
  * Option ``-g``: Generation of references only (no validation)
  * Option ``-s``: Plot differences with references only (no validation)
  * Option ``-c``: Compilation only (no run, no validation)
  * Option ``-r <nRrestarts>``: Force the simulation to be broken in several restarts.
  * Option ``-v``: Verbose
  * Option ``-h``: Help


Exit status of the script:

..
  
  * 0  validated
  * 1  validation fails
  * 2  execution fails
  * 3  compilation fails
  * 4  bad option


Examples:

..
  
  .. code-block:: bash
  
    ./validation.py -v
  
  Compiles and validates all cases in verbose mode.
  
  .. code-block:: bash
  
    ./validation.py -v -b tst1d_00_em_propagation.py 
  
  Validates only the benchmark ``tst1d_00_em_propagation.py``.
  
  .. code-block:: bash
  
    ./validation.py -v -b tst1d_00_em_propagation.py -g
  
  Generates the reference file for the benchmark ``tst1d_00_em_propagation.py``.
  
  .. code-block:: bash
  
    ./validation.py -v -b tst1d_00_em_propagation.py -s
  
  Runs the benchmark ``tst1d_00_em_propagation.py``, and plots the differences with the reference file.



.. rubric:: What does ``validation.py`` actually do?

It creates a new ``validation/workdirs`` directory (that may be freely deleted later).

It compiles the code:

..

  If the "workdirs" directory lacks a smilei binary, or it is too old,
  then the "workdirs" is backed up, and a new compilation occurs.
  The compilation output is logged in ``compilation_output``.
  If compiling errors occur, ``compilation_errors`` is created and the script exits with status 3.

It runs each benchmark:

..

  If the directory ``wd_<benchmark>/<o>/<m>`` does not exist then:
  
  * it is created.
  * ``smilei`` is executed in that directory for the requested benchmark.
  * if execution fails, the script exits with status 2.

It analyses the results (for each requested benchmark) using the ``validate_*`` script:

..

  * If requested to compare to previous references (default option), the analysis
    is compared to the reference data.
  * If requested to generate references (option ``-g``), the analysis is stored
    as reference data.
  * If requested to show differences to previous references (option ``-s``),
    the analysis is plotted vs. the reference data.


.. rubric:: How should I make the ``validate_*`` script?

The ``validate_*`` script should load the simulation results using whatever means suits
the benchmark the best. In many cases, the :doc:`happi <post-processing>` module is
employed to extract diagnostics results.

Any *python* instructions may be used to process the simulation results. Once the data
has been crunched into a meaningful value, string, or array, then it must be passed to the
following predefined function:

.. py:method:: Validate( description, data, epsilon )

  * ``description``: a string describing the data
  * ``data``: a float, a *numpy* float array, or any other python data
  * ``epsilon`` (optional): acceptable difference between data and reference

The ``data`` passed to this function constitutes the *analysis* that is compared to previous
reference files. It is the same analysis that is used to generate those reference files
in the first place.

