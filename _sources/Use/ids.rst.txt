Identification of tracked particles
-----------------------------------

:ref:`Tracked particles <DiagTrackParticles>` require an identification number (ID) in
order to be recognized after they move around. In Smilei, the particles ID are not simply
the integers from 0 to *N-1*, where *N* is the total number of particles in the
simulation. Instead, a more subtle approach is taken.

If all numbers from 0 to *N-1* were used, then processors would have to communicate
together each time a new particle is created to avoid duplicates. That would be too
costly. We choose to avoid unnecessary communications, meaning that processors manage
particle IDs independently from each other. This is reflected in the structure of the
output files for the :ref:`tracked particles diagnostic <DiagTrackParticles>`. These
files, named ``TrackParticlesDisordered_***.h5``, contain arrays where each proc *owns*
a contiguous block, corresponding to the amount of particles it needs to write::

  |------- Proc 0 ------|----------- Proc 1 ------------|--- Proc 2 ---|-- .......

These blocs have distinct sizes in general, and contain particles that are not sorted by
ID, as they move sometimes from one processor to another.

However, particles keep, in their ID, the number of the processor in which they were
created. More precisely, the ID of a particle is a ``uint64``, a positive integer whose
binary representation has 64 bits. To illustrate the content of these 64 bits let us
replace zeros and ones by X, Y or Z::
  
  XXXXXXXXYYYYYYYYYYYYYYYYYYYYYYYYZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZ

There are 64 bits available. The first 8 (X) are not parsed by the code, but may be set
by users for custom purposes. The next 24 (Y) represent the processor number. For instance,
for processor 0, all Y will equal 0; for processor 1, only the last Y will be 1. The last
32 bits (Z) indicate the particle number. This number is not unique among processors: for
example, the first particle of each proc always has number 0. The combination of these
last two numbers (YYYYY... and ZZZZZZ....) ensures a unique ID across the whole simulation.
Clearly, the IDs are not represented by a contiguous list of integers from 0 to *N-1*.
