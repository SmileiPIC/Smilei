
Physics modules
===============

The basic :doc:`algorithms` do not reproduce all the physics of a plasma.
For instance, the typical cell size is too coarse to model :doc:`collisions <collisions>`
so that a specific module is necessary. Similarly, atomic physics and
quantum processes require physics modules such as :doc:`ionization` and
:doc:`radiation_loss`.


.. toctree::
   :maxdepth: 1
    
   collisions
   ionization
   radiation_loss
   multiphoton_Breit_Wheeler