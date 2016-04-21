
#ifndef DIAGSVECTORPATCH_H
#define DIAGSVECTORPATCH_H

#include <hdf5.h>

class VectorPatch;
class Params;
class SmileiMPI;

class DiagsVectorPatch
{
public:
    static void initDumpFields    ( VectorPatch& vecPatches, Params& params, int timestep );
    static void finalizeDumpFields( VectorPatch& vecPatches, Params& params, int timestep );
    
    
    static void definePatchDiagsMaster( VectorPatch& vecPatches, hid_t globalFile, hid_t globalFileAvg );
    static void updatePatchFieldDump  ( VectorPatch& vecPatches, Params& params );


};

#endif
