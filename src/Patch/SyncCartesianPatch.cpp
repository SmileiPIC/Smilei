
#include "SyncCartesianPatch.h"

#include <vector>

#include "VectorPatch.h"
#include "Params.h"
#include "SmileiMPI.h"

using namespace std;

void SyncCartesianPatch::patchedToCartesian( VectorPatch& vecPatches, Patch* patch, Params &params, SmileiMPI* smpi, Timers &timers, int itime )
{
    SyncCartesianPatch::sync( vecPatches(0)->EMfields->Ex_, patch->EMfields->Ex_, params, smpi );

}

void SyncCartesianPatch::sync( Field* inField, Field* outField, Params &params, SmileiMPI* smpi )
{
}

