 
#include "DiagsVectorPatch.h"
#include "VectorPatch.h"
#include "Params.h"
#include "SmileiMPI.h"

#include <cstring>
#include <vector>

using namespace std;


void DiagsVectorPatch::initDumpFields(VectorPatch& vecPatches, Params& params, int timestep)
{
    vecPatches(0)->sio->createFiles(params, vecPatches(0));
    for (unsigned int ipatch=0 ; ipatch<vecPatches.size() ; ipatch++) {
        vecPatches(ipatch)->sio->setFiles( vecPatches(0)->sio->global_file_id_, vecPatches(0)->sio->global_file_id_avg );
    }
}

void DiagsVectorPatch::finalizeDumpFields(VectorPatch& vecPatches, Params& params, int timestep)
{
    for (unsigned int ipatch=1 ; ipatch<vecPatches.size() ; ipatch++) {
        vecPatches(ipatch)->sio->setFiles( 0, 0 );
    }

}

//////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////


void DiagsVectorPatch::definePatchDiagsMaster(VectorPatch& vecPatches, hid_t globalFile, hid_t globalFileAvg)
{
    for (unsigned int ipatch=0 ; ipatch<vecPatches.size() ; ipatch++) {
        vecPatches(ipatch)->sio->setFiles( globalFile, globalFileAvg );
    }

}


void DiagsVectorPatch::updatePatchFieldDump( VectorPatch& vecPatches, Params& params )
{
    for (unsigned int ipatch=0 ; ipatch<vecPatches.size() ; ipatch++) {
        if ( vecPatches(ipatch)->Pcoordinates[0]!=params.number_of_patches[0]-1 )
            vecPatches(ipatch)->sio->updatePattern( params, vecPatches(ipatch) );
    }
    
}
