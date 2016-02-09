
#include "DiagsVectorPatch.h"
#include "VectorPatch.h"
#include "Params.h"
#include "SmileiMPI.h"

#include <cstring>
#include <vector>

using namespace std;

void DiagsVectorPatch::computeGlobalDiags( VectorPatch& vecPatches, int timestep)
{
    
    DiagsVectorPatch::computeScalarsDiags  ( vecPatches, timestep );
    DiagsVectorPatch::computePhaseSpace    ( vecPatches );
    DiagsVectorPatch::computeParticlesDiags( vecPatches, timestep );
    
}

void DiagsVectorPatch::computeScalarsDiags( VectorPatch& vecPatches, int timestep)
{
    //cout << "In Global Compute Scalar Diags " << vecPatches(0)->Diags->scalars.every << " \t timestep = " << timestep << endl;
    int scalars_every( vecPatches(0)->Diags->scalars.every );
    if (timestep % scalars_every != 0) return;

    //cout << "In Global Compute Scalar Daigs\n";

    //std::vector<std::pair<std::string,double> > out_list;
    //std::vector<std::string> out_key;
    //std::vector<double>      out_value;
    //std::vector<unsigned int> out_width;
    //std::vector<std::pair<std::string,double> >::iterator itDiagScalar;


    int nDiags( vecPatches(0)->Diags->scalars.out_value.size() );
    // Initialize scalars iterator on 1st diag
    for (unsigned int ipatch=0 ; ipatch<vecPatches.size() ; ipatch++) {
        vecPatches(ipatch)->Diags->scalars.itDiagScalarName  =  vecPatches(ipatch)->Diags->scalars.out_key.begin();
        vecPatches(ipatch)->Diags->scalars.itDiagScalarValue =  vecPatches(ipatch)->Diags->scalars.out_value.begin();
    }


    for (int idiags = 0 ; idiags<nDiags ; idiags++) {
        string diagName( *vecPatches(0)->Diags->scalars.itDiagScalarName );

        if ( ( diagName.find("Min") == std::string::npos ) && ( diagName.find("Max") == std::string::npos ) ) {
            double sum(0.);
            for (unsigned int ipatch=0 ; ipatch<vecPatches.size() ; ipatch++) {
                sum += *vecPatches(ipatch)->Diags->scalars.itDiagScalarValue;
                if (ipatch) {
                    vecPatches(ipatch)->Diags->scalars.itDiagScalarName++;
                    vecPatches(ipatch)->Diags->scalars.itDiagScalarValue++;
                }
            }
            *vecPatches(0)->Diags->scalars.itDiagScalarValue = sum;
            vecPatches(0)->Diags->scalars.itDiagScalarName++;
            vecPatches(0)->Diags->scalars.itDiagScalarValue++;
        }
        else if ( diagName.find("MinCell") != std::string::npos ) {
            vector<double>::iterator iterVal    = vecPatches(0)->Diags->scalars.itDiagScalarValue-1;
            vector<double>::iterator iterValRef = vecPatches(0)->Diags->scalars.itDiagScalarValue-1;
            double min( *iterValRef );

            for (unsigned int ipatch=1 ; ipatch<vecPatches.size() ; ipatch++) {
                if ( *vecPatches(ipatch)->Diags->scalars.itDiagScalarValue < min ) {
                    min = *vecPatches(ipatch)->Diags->scalars.itDiagScalarValue;
                    iterVal = vecPatches(ipatch)->Diags->scalars.itDiagScalarValue-1;
                }
                if (ipatch) {
                    vecPatches(ipatch)->Diags->scalars.itDiagScalarName++;
                    vecPatches(ipatch)->Diags->scalars.itDiagScalarValue++;
                }
            }
            *vecPatches(0)->Diags->scalars.itDiagScalarValue = min;
            iterValRef = iterVal;

            vecPatches(0)->Diags->scalars.itDiagScalarName++;            
            vecPatches(0)->Diags->scalars.itDiagScalarValue++;            
        }
        else if ( diagName.find("MaxCell") != std::string::npos ) {
            vector<double>::iterator iterVal    = vecPatches(0)->Diags->scalars.itDiagScalarValue-1;
            vector<double>::iterator iterValRef = vecPatches(0)->Diags->scalars.itDiagScalarValue-1;
            double max( *iterValRef );

            for (unsigned int ipatch=1 ; ipatch<vecPatches.size() ; ipatch++) {
                if ( *vecPatches(ipatch)->Diags->scalars.itDiagScalarValue > max ) {
                    max = *vecPatches(ipatch)->Diags->scalars.itDiagScalarValue;
                    iterVal = vecPatches(ipatch)->Diags->scalars.itDiagScalarValue-1;
                }
                if (ipatch) {
                    vecPatches(ipatch)->Diags->scalars.itDiagScalarName++;
                    vecPatches(ipatch)->Diags->scalars.itDiagScalarValue++;
                }
            }
            *vecPatches(0)->Diags->scalars.itDiagScalarValue = max;
            iterValRef = iterVal;

            vecPatches(0)->Diags->scalars.itDiagScalarName++;            
            vecPatches(0)->Diags->scalars.itDiagScalarValue++;            
        }

        // Go to next diag
    }

    // After MPI sync
    //vecPatches(0)->Diags->scalars.write(timestep);

}

void DiagsVectorPatch::computePhaseSpace( VectorPatch& vecPatches )
{
    // A définir : DiagPhaseSpace::itDiagPhase

    int nDiags( vecPatches(0)->Diags->phases.vecDiagPhaseToRun.size() );

    // Initialize scalars iterator on 1st diag
    for (unsigned int ipatch=0 ; ipatch<vecPatches.size() ; ipatch++)
        vecPatches(ipatch)->Diags->phases.itDiagPhase =  vecPatches(ipatch)->Diags->phases.vecDiagPhaseToRun.begin();
    
    for (int idiags = 0 ; idiags<nDiags ; idiags++) {
        vector<unsigned int> diagSize = (*vecPatches(0)->Diags->phases.itDiagPhase)->my_data.dims_;
        for (unsigned int ipatch=1 ; ipatch<vecPatches.size() ; ipatch++) {
            for (int i=0 ; i<diagSize[0] ; i++)
                for (int j=0 ; j<diagSize[1] ; j++)
                    (*vecPatches(0)->Diags->phases.itDiagPhase)->my_data(i,j) += (*vecPatches(ipatch)->Diags->phases.itDiagPhase)->my_data(i,j);
            vecPatches(ipatch)->Diags->phases.itDiagPhase++;
        } // for ipatch
        vecPatches(0)->Diags->phases.itDiagPhase++;

    } // for idiags

    for (unsigned int ipatch=1 ; ipatch<vecPatches.size() ; ipatch++)
        vecPatches(ipatch)->Diags->phases.vecDiagPhaseToRun.clear();

}


void DiagsVectorPatch::computeParticlesDiags( VectorPatch& vecPatches, int timestep)
{
    int nDiags( vecPatches(0)->Diags->vecDiagnosticParticles.size() );

    for (int idiags = 0 ; idiags<nDiags ; idiags++) {
        if (timestep % vecPatches(0)->Diags->vecDiagnosticParticles[idiags]->every != vecPatches(0)->Diags->vecDiagnosticParticles[idiags]->time_average-1) continue;

        int output_size = vecPatches(0)->Diags->vecDiagnosticParticles[idiags]->output_size;
        for (unsigned int ipatch=1 ; ipatch<vecPatches.size() ; ipatch++) {
            for (int i=0 ; i<output_size ; i++)
                vecPatches(0)->Diags->vecDiagnosticParticles[idiags]->data_sum[i] += vecPatches(ipatch)->Diags->vecDiagnosticParticles[idiags]->data_sum[i];
        } // for ipatch

    } // for idiags

    
    for (unsigned int ipatch=1 ; ipatch<vecPatches.size() ; ipatch++)
        for (unsigned int i=0; i<vecPatches(ipatch)->Diags->vecDiagnosticParticles.size(); i++)
               if (vecPatches(ipatch)->Diags->vecDiagnosticParticles[i]->time_average == 1)
                   vecPatches(ipatch)->Diags->vecDiagnosticParticles[i]->clean();

}






//////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////

void DiagsVectorPatch::initProbesDiags(VectorPatch& vecPatches, Params& params, int timestep)
{
    vecPatches(0)->Diags->probes.createFile();
    // Start at 0, cause of setFile set probesStart (locate writing point in h5 file)
    for (unsigned int ipatch=0 ; ipatch<vecPatches.size() ; ipatch++) {
        vecPatches(ipatch)->Diags->probes.setFile( vecPatches(0)->Diags->probes.fileId, vecPatches(ipatch), params );
    }
    for (unsigned int ipatch=0 ; ipatch<vecPatches.size() ; ipatch++) {
        vecPatches(ipatch)->Diags->probes.waitSetFile( params );
    }    //cout << " File created " << endl;

    for (unsigned int ipatch=0 ; ipatch<vecPatches.size() ; ipatch++) {
        //cout << "Data written for " << ipatch << endl;
        vecPatches(ipatch)->Diags->probes.writePositionIn(params);
        //cout << "End of Data written for " << ipatch << endl;
    }
}

void DiagsVectorPatch::finalizeProbesDiags(VectorPatch& vecPatches, Params& params, int timestep)
{
    for (unsigned int ipatch=1 ; ipatch<vecPatches.size() ; ipatch++) {
        vecPatches(ipatch)->Diags->probes.setFile( 0 );
    }

}

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

void DiagsVectorPatch::initTrackParticles(VectorPatch& vecPatches, Params& params, SmileiMPI* smpi)
{
    int nspecies = vecPatches(0)->vecSpecies.size();
    int idiag(0);
    for ( int ispec=0 ; ispec<nspecies ; ispec++) {
        
        // Communicate some stuff if this is a species that has to be dumped (particles have Id)
        // Need to be placed after ALL createParticles()
        if (vecPatches(0)->vecSpecies[ispec]->particles->track_every) {

            // Internal patches offset

            std::vector<int> localNbrParticles( vecPatches.size(), 0 );
            localNbrParticles[0] = vecPatches(0)->vecSpecies[ispec]->getNbrOfParticles();
            for (unsigned int ipatch=1 ; ipatch<vecPatches.size() ; ipatch++) {
                // number of particles up to ipatch (including)
                localNbrParticles[ipatch] += vecPatches(ipatch)->vecSpecies[ispec]->getNbrOfParticles() + localNbrParticles[ipatch-1];
                vecPatches(ipatch)->vecSpecies[ispec]->particles->addIdOffsets(localNbrParticles[ipatch-1]);
            }
            int locNbrParticles = localNbrParticles[vecPatches.size()-1];


            // MPI offset

            //int locNbrParticles = thisSpecies->getNbrOfParticles();
            int sz(1);
            MPI_Comm_size( MPI_COMM_WORLD, &sz );
            std::vector<int> allNbrParticles(sz);
            MPI_Allgather( &locNbrParticles, 1, MPI_INTEGER, &allNbrParticles[0], 1, MPI_INTEGER, MPI_COMM_WORLD );

            int totNbrParts(0);
            for (int irk=0 ; irk<sz ; irk++) totNbrParts += allNbrParticles[irk];
            // HDF5 file open by all patch master
            vecPatches(0)->Diags->vecDiagnosticTrackParticles[idiag]->createFile(totNbrParts,params);

            // Set HDF5 context for other patches
            for (unsigned int ipatch=1 ; ipatch<vecPatches.size() ; ipatch++)
                vecPatches(ipatch)->Diags->vecDiagnosticTrackParticles[idiag]->setGlobalNbrParticles(totNbrParts);

            idiag++; // Considered DiagnosticTrackParticles ordered as Species

            int nParticles(0);

            nParticles =  allNbrParticles[0];
            for (int irk=1 ; irk<sz ; irk++){
                allNbrParticles[irk] += nParticles;
                nParticles = allNbrParticles[irk];
            }
            for (int irk=sz-1 ; irk>0 ; irk--){
                allNbrParticles[irk] = allNbrParticles[irk-1];
            }
            allNbrParticles[0] = 0;

            int offset(0);
            MPI_Scatter(&allNbrParticles[0], 1 , MPI_INTEGER, &offset, 1, MPI_INTEGER, 0, MPI_COMM_WORLD );
            
            for (unsigned int ipatch=0 ; ipatch<vecPatches.size() ; ipatch++)
                vecPatches(ipatch)->vecSpecies[ispec]->particles->addIdOffsets(offset);

        } // End if track_every

        // Count total number of track particles (need to define HDF5 context)
        
        //MPI_Bcast();

    } // End for ispec

} // End initTrackParticles



void DiagsVectorPatch::initCollisions(VectorPatch& vecPatches, Params& params, SmileiMPI* smpi)
{
    int index;
    // For each collision
    for (unsigned int icoll=0 ; icoll<vecPatches(0)->vecCollisions.size(); icoll++) {
        // The master MPI creates the debug file
        if( smpi->isMaster() ) vecPatches(0)->vecCollisions[icoll]->createFile();
        // All patch masters create arrays in the database for ionization
        index = vecPatches(0)->vecCollisions[icoll]->Ionization->createDatabase(params.wavelength_SI);
        // All patches are assigned the correct arrays in the database
        for (unsigned int ipatch=0 ; ipatch<vecPatches.size() ; ipatch++)
            vecPatches(0)->vecCollisions[icoll]->Ionization->assignDatabase(index);
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

void DiagsVectorPatch::definePatchDiagsMaster(VectorPatch& vecPatches)
{
    int patchIdMaster(0);
    for (patchIdMaster=0 ; patchIdMaster<vecPatches.size() ; patchIdMaster++ )
        if ( vecPatches(patchIdMaster)->Diags->probes.fileId != 0 ) break;

    for (unsigned int ipatch=0 ; ipatch<vecPatches.size() ; ipatch++) {
        if ((ipatch!=patchIdMaster) && (patchIdMaster!=vecPatches.size()) ) { // patchIdMaster!=vecPatches.size() 
                vecPatches(ipatch)->Diags->probes.setFile( vecPatches(patchIdMaster)->Diags->probes.fileId );
        }
    }

    for (patchIdMaster=0 ; patchIdMaster<vecPatches.size() ; patchIdMaster++ )
        if ( vecPatches(patchIdMaster)->sio->global_file_id_ != 0 ) break;

    for (unsigned int ipatch=0 ; ipatch<vecPatches.size() ; ipatch++) {
        if ((ipatch!=patchIdMaster) && (patchIdMaster!=vecPatches.size()) ) { // patchIdMaster!=vecPatches.size() 
            vecPatches(ipatch)->sio->setFiles( vecPatches(patchIdMaster)->sio->global_file_id_, vecPatches(patchIdMaster)->sio->global_file_id_avg );
        }
    }

}

void DiagsVectorPatch::updatePatchFieldDump( VectorPatch& vecPatches, Params& params )
{
    for (unsigned int ipatch=0 ; ipatch<vecPatches.size() ; ipatch++) {
        if ( vecPatches(ipatch)->Pcoordinates[0]!=params.number_of_patches[0]-1 )
            vecPatches(ipatch)->sio->updatePattern( params, vecPatches(ipatch) );
    }
    
}
