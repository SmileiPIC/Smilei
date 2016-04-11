 
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


void DiagsVectorPatch::initCollisions(VectorPatch& vecPatches, Params& params, SmileiMPI* smpi)
{
    int index;
    // For each collision
    for (unsigned int icoll=0 ; icoll<vecPatches(0)->vecCollisions.size(); icoll++) {
        // All patch masters create arrays in the database for ionization
        index = vecPatches(0)->vecCollisions[icoll]->Ionization->createDatabase(params.referenceAngularFrequency_SI);
        // All patches are assigned the correct arrays in the database
        for (unsigned int ipatch=0 ; ipatch<vecPatches.size() ; ipatch++)
            vecPatches(ipatch)->vecCollisions[icoll]->Ionization->assignDatabase(index);
    }

} // End initCollisions

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
