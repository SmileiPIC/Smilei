
#ifndef SYNCVECTORPATCH_H
#define SYNCVECTORPATCH_H

#include <vector>

class VectorPatch;
class Params;
class SmileiMPI;
class Field;
class Timers;

class SyncVectorPatch {
public :

    static void exchangeParticles(VectorPatch& vecPatches, int ispec, Params &params, SmileiMPI* smpi, Timers &timers, int itime);
    static void finalize_and_sort_parts(VectorPatch& vecPatches, int ispec, Params &params, SmileiMPI* smpi, Timers &timers, int itime);
    static void sumRhoJ  ( VectorPatch& vecPatches, Timers &timers, int itime );
    static void sumRhoJs ( VectorPatch& vecPatches, int ispec, Timers &timers, int itime );
    static void exchangeE( VectorPatch& vecPatches );
    static void finalizeexchangeE( VectorPatch& vecPatches );
    static void exchangeB( VectorPatch& vecPatches );
    static void finalizeexchangeB( VectorPatch& vecPatches );
    static void sum      ( std::vector<Field*> fields, VectorPatch& vecPatches, Timers &timers, int itime );
    static void exchange ( std::vector<Field*> fields, VectorPatch& vecPatches );
    static void finalizeexchange( std::vector<Field*> fields, VectorPatch& vecPatches );
    static void exchange0( std::vector<Field*> fields, VectorPatch& vecPatches );
    static void finalizeexchange0( std::vector<Field*> fields, VectorPatch& vecPatches );
    static void exchange1( std::vector<Field*> fields, VectorPatch& vecPatches );
    static void finalizeexchange1( std::vector<Field*> fields, VectorPatch& vecPatches );
    static void exchange2( std::vector<Field*> fields, VectorPatch& vecPatches );
    static void finalizeexchange2( std::vector<Field*> fields, VectorPatch& vecPatches );

};

#endif
