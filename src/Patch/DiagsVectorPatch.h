
#ifndef DIAGSVECTORPATCH_H
#define DIAGSVECTORPATCH_H

class VectorPatch;

class DiagsVectorPatch
{
public:
    static void computeGlobalDiags   ( VectorPatch& vecPatches, int timestep );
    static void computeScalarsDiags  ( VectorPatch& vecPatches, int timestep );
    static void computePhaseSpace    ( VectorPatch& vecPatches );
    static void computeParticlesDiags( VectorPatch& vecPatches, int timestep );

};

#endif
