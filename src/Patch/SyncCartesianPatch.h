
#ifndef SYNCCARTESIANPATCH_H
#define SYNCCARTESIANPATCH_H

class Domain;
class VectorPatch;
class Patch;
class Params;
class SmileiMPI;
class Timers;
class Field;
class ElectroMagn;

class SyncCartesianPatch
{
public :

    static void patchedToCartesian( VectorPatch &vecPatches, Domain &domain, Params &params, SmileiMPI *smpi, Timers &timers, int itime );
    static void cartesianToPatches( Domain &domain, VectorPatch &vecPatches, Params &params, SmileiMPI *smpi, Timers &timers, int itime );
    static void sync( Field *inField, Field *outField, Params &params, SmileiMPI *smpi, Patch *inPatch, Patch *outPatch );
    static void syncBack( Field *inField, Field *outField, Params &params, SmileiMPI *smpi, Patch *inPatch, Patch *outPatch );
    
};

#endif
