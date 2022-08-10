Optimization and vectorization compilation options
-----------------------------------------------------

Modern and mature compilers can automatically and deeply optimize the code for specific hardware.
Optimization features can be activated and tuned using compiler flags.
Among possible optimization capabilities, we can mention loop refactoring (unrolling, reordoring, splitting, etc), inlining, vectorization, pipelining and more.
Each compiler has its own flags.
While some share a common syntax, others provide more fine tuning options for specific processors (like Fujitsu for the A64FX).

Some optimization flags are automatically set in the makefile;
others are specified in the *machine files* because they are specific to given compilers, hardwares or processor families.
:doc:`SIMD vectorization </Understand/vectorization>` in Smilei, for instance,
can be automatic, or can rely on OpenMP's directive ``#pragma omp simd``.


Options are passed to the compiler through the environment variable ``CXXFLAGS``.
Most machine files include appropriate flags to ensure the best performance 
on a given processor type using a given compiler (most often Intel but also GNU, LLVM, ARMCLANG and Fujitsu).

For instance, :program:`Smilei` has been tested on
Intel processors (Skylake 8168) with an Intel environment.
The following flags, located in the ``skylake`` machine file, provide a good performance:

.. code-block:: bash

  -xCOMMON-AVX512 -ip -ipo -inline-factor=1000 -D__INTEL_SKYLAKE_8168 -fno-alias

The vectorization must also be activated :ref:`in the namelist <Vectorization>`.

The following explains the meaning of the flags commonly used in our machine files.
You can use this section to design your own machine file.

.. rubric:: Intel compiler flags for optimization and vectorization

The following options are commonly used by Smilei's machine files to compile on most modern x86 processors (both Intel and AMD). They enure the best performance of the code, especially for vectorization.

* ``-O3``: this option directly integrated in the makefile tells the compiler to use agressive optimization at compilation

* ``-march=cpu-type``: Tells the compiler to generate code for processors that support certain features. See `this page <https://www.intel.com/content/www/us/en/develop/documentation/cpp-compiler-developer-guide-and-reference/top/compiler-reference/compiler-options/compiler-option-details/code-generation-options/march.html>`_ for more info.

* ``-xCPU-TYPE``: Tells the compiler which processor features it may target, including which instruction sets and optimizations it may generate. ``CPU-TYPE`` can be the instruction set or the processor family name. This option is important for vectorization. See `this page <https://www.intel.com/content/www/us/en/develop/documentation/cpp-compiler-developer-guide-and-reference/top/compiler-reference/compiler-options/compiler-option-details/code-generation-options/x-qx.html#x-qx>`_ for more info.

    * ``-xCOMMON-AVX512``: for processors that support AVX512 (Skylake, Cascadelake)
    * ``-xMIC-AVX512``: suitable for first-generation AVX512 processors (KNL family)
    * ``-xCORE-AVX2``: for processors using the AVX2 instruction set such as Intel Broadwell (3nd generation E3, E5 and E7 Xeon family) and AMD Epyc processors
    * ``-xAVX``: for the first generation E3 and E5 Xeon family
    * ``-xSKYLAKE``
    * ``-xKNL``
    * ``-xBROADWELL``
    * and more...

* ``-ip``: interprocedural optimizations for single-file compilation. This flag is important for function inline in a single C++ file. This option is not by default in the makefile but is available in many machine files.
         
* ``-ipo``: Interprocedural optimization, a step that examines function calls between files when the program is linked.  This flag must be used to compile and when linking. Compile times are very long with this flag, however depending on the application there may be appreciable performance improvements when combined with the -O* flags. This flag is not by default in the makefile and is rarely used due to long compilation time. Use this flag for production runs if you do not plan to often recompile.

* ``-inline-factor=1000``: Specifies the percentage multiplier that should be applied to all inlining options that define upper limits. 
 
* ``-fno-alias``: this is recommended only as an experiment, possibly to see whether a useful restrict keyword has been missed. It would mean that the compiler could ignore unproven aliasing hazards.
 
A detailed description of the available flags for x86 instruction sets is given on this `page <https://www.intel.com/content/www/us/en/developer/articles/technical/performance-tools-compiler-options-for-sse-generation-and-processor-specific-optimizations.html/>`_.
All Intel compiler flags are listed on this  `page <https://www.intel.com/content/www/us/en/develop/documentation/cpp-compiler-developer-guide-and-reference/top/compiler-reference/compiler-options/alphabetical-list-of-compiler-options.html/>`_. 

More examples:

* For AMD EPYC ROME processors:

.. code-block:: bash

  -march=core-avx2 -xCORE-AVX2 -ip -inline-factor=1000 -fno-alias #-ipo

* For Intel Broadwell processors:

.. code-block:: bash

  -xCORE-AVX2 -O3 -ip -inline-factor=1000 -D__INTEL_BDW_E5_2697_V4 -fno-alias

* For Intel Cascadelake processors:

.. code-block:: bash

  -march=cascadelake -xCOMMON-AVX512 -ip -inline-factor=1000 -D__INTEL_CASCADELAKE_6248 -qopt-zmm-usage=high -fno-alias #-ipo

* For Intel KNL processors:

.. code-block:: bash

  -march=knl -xMIC-AVX512 -ip -inline-factor=1000 -D__INTEL_KNL_7250 -qopt-zmm-usage=high -fno-alias #-ipo

.. rubric:: GNU compiler flags for optimization and vectorization

* ``-O3``: this option directly integrated in the makefile tells the compiler to use agressive optimization at compilation

* ``-Ofast``: Disregard strict standards compliance. -Ofast enables all -O3 optimizations. It also enables optimizations that are not valid for all standard-compliant programs.It can result in incorrect output for programs that depend on an exact implementation of IEEE or ISO rules/specifications for math functions. It may, however, yield faster code for programs that do not require the guarantees of these specifications. 
    
* ``-mtune=cpu-type``: This option specifies that GCC should tune the performance of the code as if the target were of the type specified in this option, here ``cpu-type``.For some ARM implementations better performance can be obtained by using this option. Possible common cpu types are
     
    * ``cascadelake``
    * ``skylake-avx512``
    * ``a64fx`` for A64FX Fujitsu processor
    * ``knl``
    * ``broadwell``
    * ``znver2`` for 2nd generation AMD EPYC processors
    * ``znver2`` for 3rd generation AMD EPYC processors
    
* ``-march=cpu-type``: This flag does additional tuning for specific processor types. Specifying -march=cpu-type implies -mtune=cpu-type, except where noted otherwise.
    
    * ``cascadelake``
    * ``skylake-avx512``
    * ``sve`` to generate SVE instructions (vectorization on ARM like A64FX)
    * ``armv8.2-a`` to generate armv8 instructions (like A64FX)
    * ``knl``
    * ``broadwell``
    * ``znver2`` for 2nd generation AMD EPYC processors
    * ``znver3`` for 3rd generation AMD EPYC processors
    
* ``-msve-vector-bits=512``: Specify the number of bits in an SVE vector register on ARM architecture using SVE (useful for A64FX).

* ``-ffast-math``: This option is not turned on by any -O option besides -Ofast since it can result in incorrect output for programs that depend on an exact implementation of IEEE or ISO rules/specifications for math functions. It may, however, yield faster code for programs that do not require the guarantees of these specifications. 

Find out more information:

* `Optimization flags for GNU <https://gcc.gnu.org/onlinedocs/gcc/Optimize-Options.html>`_.
* `g++ man page <https://linux.die.net/man/1/g++>`_.
* `list of architecture options <https://gcc.gnu.org/onlinedocs/gcc/Submodel-Options.html#Submodel-Options/>`_.
* `AMD quick reference guide <https://developer.amd.com/wordpress/media/2020/04/Compiler%20Options%20Quick%20Ref%20Guide%20for%20AMD%20EPYC%207xx2%20Series%20Processors.pdf>`_.

.. rubric:: LLVM compiler flags for optimization and vectorization

The LLVM compiler shares many flags with the GNU one.

* ``-O3``: this option directly integrated in the makefile tells the compiler to use agressive optimization at compilation

* ``-Ofast``: Disregard strict standards compliance. -Ofast enables all -O3 optimizations. It also enables optimizations that are not valid for all standard-compliant programs.It can result in incorrect output for programs that depend on an exact implementation of IEEE or ISO rules/specifications for math functions. It may, however, yield faster code for programs that do not require the guarantees of these specifications. 

* ``-mtune=cpu-type``: This option specifies that LLVM should tune the performance of the code as if the target were of the type specified in this option, here ``cpu-type``.For some ARM implementations better performance can be obtained by using this option. Possible common cpu types are
     
    * ``cascadelake``
    * ``skylake-avx512``
    * ``a64fx`` for A64FX Fujitsu processor
    * ``knl``
    * ``broadwell``
    * ``znver2`` for 2nd generation AMD EPYC processors
    * ``znver2`` for 3rd generation AMD EPYC processors
    
* ``-march=cpu-type``: This flag does additional tuning for specific processor types. Specifying ``-march=cpu-type`` implies ``-mtune=cpu-type``, except where noted otherwise.
    
    * ``cascadelake``
    * ``skylake-avx512``
    * ``sve`` to generate SVE instructions (vectorization on ARM like A64FX)
    * ``armv8.2-a`` to generate armv8 instructions (like A64FX)
    * ``knl``
    * ``broadwell``
    * ``znver2`` for 2nd generation AMD EPYC processors
    * ``znver3`` for 3rd generation AMD EPYC processors

* ``-ffast-math``: Enable fast-math mode. This option lets the compiler make aggressive, potentially-lossy assumptions about floating-point math.

* ``-ffinite-math-only``: Allow floating-point optimizations that assume arguments and results are not NaNs or +-Inf.

* ``-ffp-contract=off|on|fast|fast-honor-pragmas``: Specify when the compiler is permitted to form fused floating-point operations, such as fused multiply-add (FMA). Fused operations are permitted to produce more precise results than performing the same operations separately.

Some examples:

* For Intel Cascade processors:

.. code-block:: bash

  -mtune=cascadelake -march=cascadelake -ffinite-math-only -ffp-contract=fast

* For the A64FX processor:

.. code-block:: bash

  -march=armv8.2-a+sve -ffinite-math-only -fsimdmath -fopenmp-simd -ffp-contract=fast #-ffast-math

Find out more information:

* `LLVM Clang user manual <https://clang.llvm.org/docs/UsersManual.html>`_.

.. rubric:: Fujitsu compiler flags for optimization and vectorization

Fujitsu compiler is only used or the A64FX processor so far.
The compiler can work in two differents modes called Trad and Clang mode.
The Clang mode uses the Clang compilation flags.
The Trad moe is usually the default one, the Clang mode an be activated using the flag ``-Nclang``.

* ``-O3`` (bothy in Trad and clang mode): by default in the makefile

* ``-Kfast`` (in Trad mode) / ``-Ofast`` (in Clang mode)

* ``-KA64FX`` (in Trad mode)

* ``-KSVE`` (in Trad mode) / ``-march=sve`` (in Clang mode)

* ``-KSSL2`` (in Trad mode)

* ``-Kparallel`` (in Trad mode)

* ``-Kunroll`` (in Trad mode)

* ``-Ksimd=2`` (in Trad mode)

* ``-Kassume=notime_saving_compilation`` (in Trad mode)

* ``-Kocl`` (in Trad mode) / ``-ffj-ocl`` (in Clang mode)

.. rubric:: Smilei compiler flags for adaptive vectorization

Performance models are implemented in :program:`Smilei` for adaptive vectorization.
By default, a general performance model is used but some performance models can be used for specific types of processors:

* ``-D__INTEL_CASCADELAKE_6248``
* ``-D__INTEL_SKYLAKE_8168``
* ``-D__AMD_ROME_7H12``
* ``-D__INTEL_KNL_7250``: available in 3D only
* ``-D__INTEL_BDW_E5_2697_V4``: available in 3D only
* ``-D__INTEL_HSW_E5_2680_v3``: available in 3D only

These flags are used in the corresponding machine files.
