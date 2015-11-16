/*
 * SmileiIO.cpp
 *
 *  Created on: 3 juil. 2013
 */

#include "SmileiIO.h"

#include <sstream>
#include <iomanip>

#include <mpi.h>

#include "Params.h"
#include "Diagnostic.h"
#include "Patch.h"
#include "SmileiMPI.h"
#include "SimWindow.h"
#include "ElectroMagn.h"
#include "Species.h"


using namespace std;

SmileiIO::SmileiIO( Params& params, Diagnostic *diag, Patch* patch ) : 
fieldsToDump(diag->fieldsToDump),
global_file_id_(0),
global_file_id_avg(0)
{
    //
    // Create property list for collective dataset write.
    //
    write_plist = H5Pcreate(H5P_DATASET_XFER);
    H5Pset_dxpl_mpio(write_plist, H5FD_MPIO_INDEPENDENT);
}


void SmileiIO::setFiles( hid_t masterFileId, hid_t masterFileIdAvg )
{
    global_file_id_ = masterFileId;
    global_file_id_avg = masterFileIdAvg;
}

void SmileiIO::createFiles( Params& params, Patch* patch)
{
    
    // ----------------------------
    // Management of global IO file
    // ----------------------------
    MPI_Info info  = MPI_INFO_NULL;
    hid_t plist_id = H5Pcreate(H5P_FILE_ACCESS);
    H5Pset_fapl_mpio(plist_id, MPI_COMM_WORLD, info);

    // Fields.h5
    // ---------
    global_file_id_  = H5Fcreate( "Fields.h5",     H5F_ACC_TRUNC, H5P_DEFAULT, plist_id);

    // Create property list for collective dataset write: for Fields.h5
    H5::attr(global_file_id_, "res_time", params.res_time);
    H5::attr(global_file_id_, "every", patch->Diags->fieldDump_every);
    H5::attr(global_file_id_, "res_space", params.res_space);
    H5::attr(global_file_id_, "sim_length", params.sim_length);

    // Fields_avg.h5
    // -------------
    global_file_id_avg = 0;
    if  (patch->Diags->ntime_step_avg!=0) {
        global_file_id_avg = H5Fcreate( "Fields_avg.h5", H5F_ACC_TRUNC, H5P_DEFAULT, plist_id);
        
        // Create property list for collective dataset write: for Fields.h5
        H5::attr(global_file_id_avg, "res_time", params.res_time);
        H5::attr(global_file_id_avg, "every", patch->Diags->fieldDump_every);
        H5::attr(global_file_id_avg, "res_space", params.res_space);
        H5::attr(global_file_id_avg, "sim_length", params.sim_length);
    }
    
    H5Pclose(plist_id);
    
}

SmileiIO::~SmileiIO()
{
    // Management of global IO file
    if (global_file_id_ != 0)
	H5Fclose( global_file_id_ );
    // Management of global IO file
    if (global_file_id_avg != 0)
        H5Fclose( global_file_id_avg );

    H5Pclose( write_plist );
}


void SmileiIO::createTimeStepInSingleFileTime( int time, Diagnostic* diag )
{
    ostringstream name_t;
    name_t.str("");
    name_t << "/" << setfill('0') << setw(10) << time;
	
    DEBUG(10,"[hdf] GROUP _________________________________ " << name_t.str());
    hid_t group_id = H5Gcreate(global_file_id_, name_t.str().c_str(), H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
    H5Gclose(group_id);

    if (diag->ntime_step_avg!=0) {
	group_id = H5Gcreate(global_file_id_avg, name_t.str().c_str(), H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
	H5Gclose(group_id);
    }
 
}

// ---------------------------------------------------------------------------------------------------------------------
// Write all fields of all time step in the same file
// ---------------------------------------------------------------------------------------------------------------------
void SmileiIO::writeAllFieldsSingleFileTime( std::vector<Field*> * fields, int time, bool avg )
{
    // Make group name: "/0000000000", etc.
    ostringstream name_t;
    name_t.str("");
    name_t << "/" << setfill('0') << setw(10) << time;

    // Create group inside HDF5 file
    hid_t file_id;
    if( avg ) file_id = global_file_id_avg; // different file for avg fields
    else      file_id = global_file_id_;
    //hid_t group_id = H5::group(file_id, name_t.str());
    hid_t group_id = H5Gopen(file_id, name_t.str().c_str(), H5P_DEFAULT);

    int nFields = fields->size();
    int nFieldsToDump = fieldsToDump.size();
    
    if (!nFieldsToDump)
        for (int i=0; i<nFields; i++)
            writeFieldsSingleFileTime( (*fields)[i], group_id );
    else
        for (int i=0; i<nFields; i++)
            for (int j=0; j<nFieldsToDump; j++)
                if ((*fields)[i]->name==fieldsToDump[j])
                    writeFieldsSingleFileTime( (*fields)[i], group_id );
    
    H5Gclose(group_id);
    
    H5Fflush( file_id, H5F_SCOPE_GLOBAL );
    
}


void SmileiIO::initWriteTestParticles(Species* species, int ispec, int time, Params& params, SmileiMPI* smpi) {

    Particles &cuParticles = (*species->particles);
    int locNbrParticles = species->getNbrOfParticles();
    
    hsize_t nParticles = (hsize_t) smpi->globalNbrParticles(species, locNbrParticles);
    
    if ( smpi->isMaster() && true ) {
    
        ostringstream nameDump("");
        nameDump << "TestParticles_" << species->species_param.species_type  << ".h5" ;
        hid_t fid = H5Fcreate( nameDump.str().c_str(), H5F_ACC_TRUNC, H5P_DEFAULT, H5P_DEFAULT);
        
        
        hsize_t dimsPart[2] = {0, nParticles};
        hsize_t maxDimsPart[2] = {H5S_UNLIMITED, nParticles};
        
        hid_t file_space = H5Screate_simple(2, dimsPart, maxDimsPart);
        
        
        hid_t plist = H5Pcreate(H5P_DATASET_CREATE);
        H5Pset_layout(plist, H5D_CHUNKED);
        hsize_t chunk_dims[2] = {1, nParticles};
        H5Pset_chunk(plist, 2, chunk_dims);
        
        
        hid_t did;
        for (unsigned int i=0; i<cuParticles.Position.size(); i++) {
            ostringstream namePos("");
            namePos << "Position-" << i;
            did = H5Dcreate(fid, namePos.str().c_str(), H5T_NATIVE_DOUBLE, file_space, H5P_DEFAULT, plist, H5P_DEFAULT);
            H5Dclose(did);
        }
        
        for (unsigned int i=0; i<cuParticles.Momentum.size(); i++) {
            ostringstream namePos("");
            namePos << "Momentum-" << i;
            did = H5Dcreate(fid, namePos.str().c_str(), H5T_NATIVE_DOUBLE, file_space, H5P_DEFAULT, plist, H5P_DEFAULT);
            H5Dclose(did);
        }
        
        did = H5Dcreate(fid, "Weight", H5T_NATIVE_DOUBLE, file_space, H5P_DEFAULT, plist, H5P_DEFAULT);
        H5Dclose(did);
    
        did = H5Dcreate(fid, "Weight", H5T_NATIVE_DOUBLE, file_space, H5P_DEFAULT, plist, H5P_DEFAULT);
        H5Dclose(did);
        
        did = H5Dcreate(fid, "Charge", H5T_NATIVE_SHORT, file_space, H5P_DEFAULT, plist, H5P_DEFAULT);
        H5Dclose(did);
        
        
        if (cuParticles.isTestParticles) {
            did = H5Dcreate(fid, "Id", H5T_NATIVE_UINT, file_space, H5P_DEFAULT, plist, H5P_DEFAULT);
            H5Dclose(did);
            
            
            if (cuParticles.isTestParticles) {
                did = H5Dcreate(fid, "Id", H5T_NATIVE_SHORT, file_space, H5P_DEFAULT, plist, H5P_DEFAULT);
                H5Dclose(did);
            }
            
            H5Pclose(plist);
            H5Sclose(file_space);
            
            
            H5Fclose( fid );
        }
    }
}

void SmileiIO::writeTestParticles(Species* species, int ispec, int time, Params& params, SmileiMPI* smpi) {
    
    // Master gathers the number of test particles
    Particles &cuParticles = (*species->particles);
    int locNbrParticles = species->getNbrOfParticles();
    int* allNbrParticles = new int[smpi->smilei_sz];
    MPI_Gather( &locNbrParticles, 1, MPI_INTEGER, allNbrParticles, 1, MPI_INTEGER, 0, smpi->SMILEI_COMM_WORLD );
    
    // Create temporary Particles object
    Particles testParticles;
    int nTestParticles(0);
    testParticles.initialize( 0, params, ispec);
    
    // Master gathers the test particles from all procs in testParticles
    if ( smpi->isMaster() ) {
        MPI_Status status;
        cuParticles.cp_particles( 0, allNbrParticles[0], testParticles , 0);
        nTestParticles+=allNbrParticles[0];
        for (int irk=1 ; irk<smpi->getSize() ; irk++) {
            if (allNbrParticles[irk]!=0) {
                Particles partVectorRecv;
                partVectorRecv.initialize( allNbrParticles[irk], params, ispec );
                MPI_Datatype typePartRecv = smpi->createMPIparticles( &partVectorRecv, params.nDim_particle + 3 + 1 + 1 );
                MPI_Recv( &(partVectorRecv.position(0,0)), 1, typePartRecv,  irk, irk, smpi->SMILEI_COMM_WORLD, &status );
                MPI_Type_free( &typePartRecv );
                // cp in testParticles
                partVectorRecv.cp_particles( 0, allNbrParticles[irk], testParticles ,nTestParticles );
                nTestParticles+=allNbrParticles[irk];
            }
        }
        
    }
    // Other procs send their test particles to Master
    else if ( locNbrParticles ) {
        MPI_Datatype typePartSend = smpi->createMPIparticles( &cuParticles, params.nDim_particle + 3 + 1 + 1 );
        MPI_Send( &(cuParticles.position(0,0)), 1, typePartSend,  0, smpi->getRank(), smpi->SMILEI_COMM_WORLD );
        MPI_Type_free( &typePartSend );
    }
    delete [] allNbrParticles;
    
    // Master writes the particles
    if ( smpi->isMaster()  ) {
        
        // Sort test particles before dump
        testParticles.sortById();
        
        ostringstream nameDump("");
        nameDump << "TestParticles_" << species->species_param.species_type  << ".h5" ;
        hid_t fid = H5Fopen( nameDump.str().c_str(), H5F_ACC_RDWR, H5P_DEFAULT);                        
        
        ostringstream attr("");
        for (unsigned int idim=0 ; idim<params.nDim_particle ; idim++) {
            attr.str("");
            attr << "Position-" << idim;
            appendTestParticles0( fid, attr.str(), testParticles.position(idim), nTestParticles, H5T_NATIVE_DOUBLE );
        }
        for (unsigned int idim=0 ; idim<3 ; idim++) {
            attr.str("");
            attr << "Momentum-" << idim;
            appendTestParticles0( fid, attr.str(), testParticles.momentum(idim), nTestParticles, H5T_NATIVE_DOUBLE );
        }
        appendTestParticles0( fid, "Weight", testParticles.weight(), nTestParticles, H5T_NATIVE_DOUBLE );
        appendTestParticles0( fid, "Charge", testParticles.charge(), nTestParticles, H5T_NATIVE_SHORT );
        if (species->particles->isTestParticles)
            appendTestParticles0( fid, "Id", testParticles.id(), nTestParticles, H5T_NATIVE_UINT );
        
        H5Fclose( fid );
        
    }
    smpi->barrier();

}


template <class T>
void SmileiIO::appendTestParticles0( hid_t fid, string name, std::vector<T> property, int nParticles, hid_t type) {

    hsize_t count[2] = {1, (hsize_t)nParticles};
    hid_t partMemSpace = H5Screate_simple(2, count, NULL);
    
    hid_t did = H5Dopen( fid, name.c_str(), H5P_DEFAULT );
    
    hid_t file_space = H5Dget_space(did);
    hsize_t dimsO[2];
    H5Sget_simple_extent_dims(file_space, dimsO, NULL);
    H5Sclose(file_space);
    hsize_t dims[2];
    dims[0] = dimsO[0]+1;
    dims[1] = dimsO[1];
    H5Dset_extent(did, dims);
    file_space = H5Dget_space(did);
    hsize_t start[2];
    start[0] = dimsO[0];
    start[1] = 0;
    H5Sselect_hyperslab(file_space, H5S_SELECT_SET, start, NULL, count, NULL);
    H5Dwrite( did, type, partMemSpace, file_space, H5P_DEFAULT, &(property[0]) );
    H5Sclose(file_space);
    
    
    H5Sclose(partMemSpace);
    H5Dclose(did);
    H5Fflush( fid, H5F_SCOPE_GLOBAL );

}
