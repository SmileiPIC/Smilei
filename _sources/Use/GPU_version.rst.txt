GPU version of SMILEI
-------


This page contains the links of this documentation to compile and run SMILEI on GPU clusters, as well as the list of the features that are supported or in development

* :doc:`Introduction</Overview/highlights>` #Improved performance using GPU offloading

* :doc:`GPU offloading , supported features, guidelines and roadmap</Understand/GPU_offloading>`

* :doc:`Intallation and compilation</Use/install_linux_GPU>`

* :doc:`Compilation for GPU-accelerated nodes</Use/installation>` #compilation-for-gpu-accelerated-nodes

* :doc:`Running on GPU-equiped nodes</Use/run>` #running-on-gpu-equiped-nodes

----

Important note: 

The biggest challenge to execute SMILEI on an accelerator is the correct installation of the openmpi library. It needs to be compiled with nvc++ after configuring (ie. ./configure --options) with the appropriate options specific to your system 
