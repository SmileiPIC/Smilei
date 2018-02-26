
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
    static void sumRhoJ  ( Params& params, VectorPatch& vecPatches, Timers &timers, int itime );
    static void sumRhoJs ( Params& params, VectorPatch& vecPatches, int ispec, Timers &timers, int itime );
    static void exchangeE( Params& params, VectorPatch& vecPatches );
    static void finalizeexchangeE( Params& params, VectorPatch& vecPatches );
    static void exchangeB( Params& params, VectorPatch& vecPatches );
    static void exchangeJ( Params& params, VectorPatch& vecPatches );
    static void exchangeA( Params& params, VectorPatch& vecPatches );
    static void exchangeGradPhi( Params& params, VectorPatch& vecPatches );
    static void finalizeexchangeB( Params& params, VectorPatch& vecPatches );
    static void finalizeexchangeJ( Params& params, VectorPatch& vecPatches );
    static void finalizeexchangeA( Params& params, VectorPatch& vecPatches );
    static void finalizeexchangeGradPhi( Params& params, VectorPatch& vecPatches );
    

    static void sumRhoJ  ( Params& params, VectorPatch& vecPatches, int imode, Timers &timers, int itime );
    static void exchangeB( Params& params, VectorPatch& vecPatches, int imode );
    static void finalizeexchangeB( Params& params, VectorPatch& vecPatches, int imode );

    static void sum       ( std::vector<Field*> fields, VectorPatch& vecPatches, Timers &timers, int itime );
    static void sumComplex( std::vector<Field*> fields, VectorPatch& vecPatches, Timers &timers, int itime );
    static void new_sum      ( std::vector<Field*>& fields, VectorPatch& vecPatches, Timers &timers, int itime );
    static void exchange ( std::vector<Field*> fields, VectorPatch& vecPatches );
    static void finalizeexchange( std::vector<Field*> fields, VectorPatch& vecPatches );
    static void exchangeComplex ( std::vector<Field*> fields, VectorPatch& vecPatches );
    static void finalizeexchangeComplex( std::vector<Field*> fields, VectorPatch& vecPatches );

    static void exchange_per_direction ( std::vector<Field*> fields, VectorPatch& vecPatches );

    static void exchange0( std::vector<Field*> fields, VectorPatch& vecPatches );
    static void new_exchange0( std::vector<Field*>& fields, VectorPatch& vecPatches );
    static void finalizeexchange0( std::vector<Field*> fields, VectorPatch& vecPatches );
    static void new_finalizeexchange0( std::vector<Field*>& fields, VectorPatch& vecPatches );
    static void exchange1( std::vector<Field*> fields, VectorPatch& vecPatches );
    static void new_exchange1( std::vector<Field*>& fields, VectorPatch& vecPatches );
    static void finalizeexchange1( std::vector<Field*> fields, VectorPatch& vecPatches );
    static void new_finalizeexchange1( std::vector<Field*>& fields, VectorPatch& vecPatches );
    static void exchange2( std::vector<Field*> fields, VectorPatch& vecPatches );
    static void new_exchange2( std::vector<Field*> fields, VectorPatch& vecPatches );
    static void finalizeexchange2( std::vector<Field*> fields, VectorPatch& vecPatches );
    static void new_finalizeexchange2( std::vector<Field*> fields, VectorPatch& vecPatches );

};

#endif
