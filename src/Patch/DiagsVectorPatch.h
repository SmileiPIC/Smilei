
#ifndef DIAGSVECTORPATCH_H
#define DIAGSVECTORPATCH_H

#include <hdf5.h>

class VectorPatch;
class Params;
class SmileiMPI;

class DiagsVectorPatch
{
public:
    static void computeGlobalDiags   ( VectorPatch& vecPatches, int timestep );
    static void computeScalarsDiags  ( VectorPatch& vecPatches, int timestep );
    static void computePhaseSpace    ( VectorPatch& vecPatches );
    static void computeParticlesDiags( VectorPatch& vecPatches, int timestep );
    
    static void initProbesDiags    ( VectorPatch& vecPatches, Params& params, int timestep );
    static void finalizeProbesDiags( VectorPatch& vecPatches, Params& params, int timestep );
    
    static void initDumpFields    ( VectorPatch& vecPatches, Params& params, int timestep );
    static void finalizeDumpFields( VectorPatch& vecPatches, Params& params, int timestep );
    
    static void initTrackParticles( VectorPatch& vecPatches, Params& params, SmileiMPI* smpi );
    
    static void initCollisions( VectorPatch& vecPatches, Params& params, SmileiMPI* smpi );
    
    static void definePatchDiagsMaster( VectorPatch& vecPatches, hid_t globalFile, hid_t globalFileAvg );
    static void definePatchDiagsMaster( VectorPatch& vecPatches );
    static void updatePatchFieldDump  ( VectorPatch& vecPatches, Params& params );


};

#endif
