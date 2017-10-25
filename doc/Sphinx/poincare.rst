

Hands-on environment
--------------------

I. Station
^^^^^^^^^^
* login : ``training``
* password : ``s@clay``

1. Download from github

Note : a SSh key had been added to a github account (Settings/SSH and GPG keys)

Can be download directly on : https://github.com/SmileiPIC/Smilei/archive/master.zip

.. code-block:: bash

   [training@mdlspc140:~/] $ git clone git@github.com:SmileiPIC/Smilei.git
   Clonage dans 'Smilei'...
   ...
   Résolution des deltas: 100% (21592/21592), fait.
   Vérification de la connectivité... fait.
   [training@mdlspc140:~/] $ cd Smilei

2. Compile the documentation or obtain online
   
.. code-block:: bash

   [training@mdlspc140:~/Smilei] $ make doc
   Compiling sphinx documentation in doc/html/Sphinx/html
   make[1] : on entre dans le répertoire « /home/training/Smilei/doc/Sphinx »
   sphinx-build -b html -d ../html/Sphinx/doctrees   . ../html/Sphinx/html
   Running Sphinx v1.3.6
   making output directory...
   ...
   finished...
   [training@mdlspc140:~/Smilei] $  firefox doc/html/Sphinx/html/index.html &                 


II. Poincare
^^^^^^^^^^^^

1. Upload Smilei on Poincare

.. code-block:: bash

   [training@mdlspc140:~/Smilei] $ cd
   [training@mdlspc140:~/] $ scp -r Smilei poincare:~/
   [training@mdlspc140:~/] $ ssh poincare -X
   ...             
   [training01@poincareint01:~/] $
              

2. Resources description

https://groupes.renater.fr/wiki/poincare/public/description_de_poincare

Node :
   
* Compute : 2 Xeon Sandy Bridge of 8 cores
* Memory : 32 Go

Software environment :
       
* compiler : Intel
* MPI : IntelMPI (``MPI_THREAD_MULTIPLE``)
* HDF5 : compiled using IntelMPI
* GNU : C++11 compatible
* Anaconda : rich Python distribution for post-processing 

Yet loaded in your environment :

.. code-block:: bash

   [training01@poincareint01:~/] $ module list
   Currently Loaded Modulefiles:
      1) intel/15.0.0                    3) hdf5/1.8.16_intel_intelmpi_mt   5) gnu/4.7.2
      2) intelmpi/5.0.1                  4) python/anaconda-2.1.0
          
Smilei
^^^^^^

1. Compile the code

.. code-block:: bash

   [training01@poincareint01:~/] $ cd Smilei
   [training01@poincareint01:~/Smilei] $ make -j4
   Checking dependencies for src/Tools/Timer.cpp
   Checking dependencies for src/Tools/tabulatedFunctions.cpp
   ...
   Compiling src/Checkpoint/Checkpoint.cpp
   Compiling src/Collisions/CollisionalIonization.cpp
   ...
   Compiling src/Tools/Timer.cpp
   Linking smilei
   Compiling src/Smilei.cpp for test mode
   Linking smilei_test for test mode
   [training01@poincareint01~/Smilei] $ ls smilei smilei_test
   smilei  smilei_test

2.  Test smilei
      
.. code-block:: bash

   [training01@poincare026-adm:~/Smilei] cp ... test.py
   [training01@poincare026-adm:~/Smilei] ./smilei_test test.py
   ...

3.  Execute smilei

Set minimal OpenMP runtime environment :

.. code-block:: bash

   [training01@poincareint01:~/Smilei] $ cat scripts/set_omp_env.sh
   #!/bin/bash

   export OMP_NUM_THREADS=$1
   export OMP_SCHEDULE=dynamic
   export OMP_PROC_BIND=true

   [training01@poincareint01:~/Smilei] $ . scripts/set_omp_env.sh 8

Single node :
   
.. code-block:: bash

   [training01@poincareint01:~/Smilei] $ llinteractive 1 clallmds+ 2
   [training17@poincare026-adm:~/Smilei] $ llq -j $LOADL_JOB_NAME
   Id                       Owner      Submitted   ST PRI Class        Running On 
   ------------------------ ---------- ----------- -- --- ------------ -----------
   poincareint02-adm.25621- training17 10/24 13:09 R  50  clallmds     poincare026-adm
   
   1 job step(s) in query, 0 waiting, 0 pending, 1 running, 0 held, 0 preempted
   [training01@poincare026-adm:~/Smilei] mpirun -np 2 ./smilei test.py
   ...             
   [training01@poincare026-adm:~/Smilei] $ exit
   logout
   Connection to poincare026-adm.maisondelasimulation.fr closed.

Multi nodes :
   
.. code-block:: bash

   [training01@poincareint02:~/Smilei] $  llinteractive 2 clallmds+ 2
   [training01@poincare026-adm:~/Smilei] $ $ llnodes.py $LOADL_JOB_NAME 
   poincareint02-adm.maisondelasimulation.fr.25622  :  2 -  poincare[026-027]
   [training01@poincare026-adm:~/Smilei] mpirun -np 4 -ppn 2 -print-rank-map ./smilei test.py
   ...             
   [training01@poincare026-adm:~/Smilei] $ exit
   logout
   Connection to poincare026-adm.maisondelasimulation.fr closed.

4. Post-processing

.. code-block:: bash

    [training01@poincareint02:~/Smilei] $ make install_python
    [training01@poincareint02:~/Smilei] $ ipython
    In [1]: import happi
