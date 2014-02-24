/*
 * SmileiIO.cpp
 *
 *  Created on: 3 juil. 2013
 *      Author: jderouil
 */

#include "SmileiIO.h"

#include "PicParams.h"
#include "SmileiMPI.h"
#include "ElectroMagn.h"
#include "Species.h"

#include <sstream>
#include <iomanip>

using namespace std;

SmileiIO::SmileiIO( PicParams* params, SmileiMPI* smpi )
{
    ostringstream name("");
    name << "particles-" << setfill('0') << setw(4) << smpi->getRank() << ".h5" ;

    hid_t attribute_id;

    // Create 1 file containing 1 dataset per Species
    partFile_id = H5Fcreate( name.str().c_str(), H5F_ACC_TRUNC, H5P_DEFAULT, H5P_DEFAULT);

    particleSize = params->nDim_particle + 3 + 1;
    hsize_t dims[2] = {0, particleSize};
    hsize_t max_dims[2] = {H5S_UNLIMITED, particleSize};
    hid_t file_space = H5Screate_simple(2, dims, max_dims);

    hid_t plist = H5Pcreate(H5P_DATASET_CREATE);
    H5Pset_layout(plist, H5D_CHUNKED);
    hsize_t chunk_dims[2] = {1, particleSize};
    H5Pset_chunk(plist, 2, chunk_dims);

    nDatasetSpecies = params->n_species;
    partDataset_id = new hid_t[nDatasetSpecies];
    for (unsigned int ispec=0 ; ispec<nDatasetSpecies ; ispec++) {
        ostringstream speciesName("");
        speciesName << params->species_param[ispec].species_type;
        partDataset_id[ispec] = H5Dcreate(partFile_id, speciesName.str().c_str(), H5T_NATIVE_FLOAT, file_space, H5P_DEFAULT, plist, H5P_DEFAULT);

        hsize_t scalaire = 1;
        hid_t tmp_space = H5Screate_simple(1, &scalaire, NULL);

        attribute_id = H5Acreate2 (partDataset_id[ispec], "Mass", H5T_IEEE_F64BE, tmp_space,
                                   H5P_DEFAULT, H5P_DEFAULT);
        H5Awrite(attribute_id, H5T_NATIVE_DOUBLE, &params->species_param[ispec].mass);
        H5Aclose(attribute_id);

        attribute_id = H5Acreate2 (partDataset_id[ispec], "Charge", H5T_IEEE_F64BE, tmp_space,
                                   H5P_DEFAULT, H5P_DEFAULT);
        H5Awrite(attribute_id, H5T_NATIVE_DOUBLE, &params->species_param[ispec].charge);
        H5Aclose(attribute_id);

        H5Sclose(tmp_space);
    }

    H5Pclose(plist);
    H5Sclose(file_space);

    dims[0] = 1;
    dims[1] = particleSize;
    partMemSpace = H5Screate_simple(2, dims, NULL);

    // Management of global IO file
    MPI_Info info  = MPI_INFO_NULL;
    hid_t plist_id = H5Pcreate(H5P_FILE_ACCESS);
    H5Pset_fapl_mpio(plist_id, MPI_COMM_WORLD, info);
    global_file_id_ = H5Fcreate( "Fields.h5", H5F_ACC_TRUNC, H5P_DEFAULT, plist_id);
    H5Pclose(plist_id);

    //
    // Create property list for collective dataset write.
    //
    write_plist = H5Pcreate(H5P_DATASET_XFER);
    H5Pset_dxpl_mpio(write_plist, H5FD_MPIO_INDEPENDENT);


}

SmileiIO::~SmileiIO()
{
    // Management of global IO file
    H5Fclose( global_file_id_ );

    H5Sclose(partMemSpace);
    for ( unsigned int s=0 ; s<nDatasetSpecies ; s++ )
        H5Dclose(partDataset_id[s]);
    delete [] partDataset_id;
    H5Fclose(partFile_id);

    H5Pclose( write_plist );
}

// ---------------------------------------------------------------------------------------------------------------------
// Write all fields of all time step in the same file
// ---------------------------------------------------------------------------------------------------------------------
void SmileiIO::writeAllFieldsSingleFileTime( ElectroMagn* EMfields, int time )
{
    ostringstream name_t;
    name_t.str("");
    name_t << "/" << setfill('0') << setw(10) << time;

    DEBUG(10,"[hdf] GROUP _________________________________ " << name_t.str());
    hid_t group_id = H5Gcreate2(global_file_id_, name_t.str().c_str(), H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);

    writeFieldsSingleFileTime( EMfields->Ex_, group_id );
    writeFieldsSingleFileTime( EMfields->Ey_, group_id );
    writeFieldsSingleFileTime( EMfields->Ez_, group_id );
    writeFieldsSingleFileTime( EMfields->Bx_m, group_id );
    writeFieldsSingleFileTime( EMfields->By_m, group_id );
    writeFieldsSingleFileTime( EMfields->Bz_m, group_id );
    writeFieldsSingleFileTime( EMfields->Jx_, group_id );
    writeFieldsSingleFileTime( EMfields->Jy_, group_id );
    writeFieldsSingleFileTime( EMfields->Jz_, group_id );
    writeFieldsSingleFileTime( EMfields->rho_, group_id );

    // for all species related quantities
    for (unsigned int ispec=0; ispec<EMfields->n_species; ispec++) {
        writeFieldsSingleFileTime( EMfields->rho_s[ispec], group_id );
        writeFieldsSingleFileTime( EMfields->Jx_s[ispec],  group_id );
        writeFieldsSingleFileTime( EMfields->Jy_s[ispec],  group_id );
        writeFieldsSingleFileTime( EMfields->Jz_s[ispec],  group_id );
    }


    H5Gclose(group_id);

    H5Fflush( global_file_id_, H5F_SCOPE_GLOBAL );

}


// ---------------------------------------------------------------------------------------------------------------------
// Each MPI process writes is particles in its own file, data are overwritten at each call ( particles-MPI_Rank.h5 )
// In progress ...
// ---------------------------------------------------------------------------------------------------------------------
void SmileiIO::writePlasma( vector<Species*> vecSpecies, double time, SmileiMPI* smpi )
{

    if (smpi->isMaster()) MESSAGE("write species disabled");
    return;
    int n_species = vecSpecies.size();
    for (int ispec=0 ; ispec<n_species ; ispec++) {
        Particles* cuParticles = &(vecSpecies[ispec])->particles;
        MESSAGE(2,"write species " << ispec);

        for (unsigned int p=0; p<(vecSpecies[ispec])->getNbrOfParticles(); p++ ) {

            hid_t file_space = H5Dget_space(partDataset_id[ispec]);
            hsize_t dimsO[2];
            H5Sget_simple_extent_dims(file_space, dimsO, NULL);
            H5Sclose(file_space);
            hsize_t dims[2];
            dims[0] = dimsO[0]+1;
            dims[1] = dimsO[1];
            H5Dset_extent(partDataset_id[ispec], dims);

            file_space = H5Dget_space(partDataset_id[ispec]);
            hsize_t start[2];
            hsize_t count[2] = {1, particleSize};
            start[0] = dimsO[0];
            start[1] = 0;
            H5Sselect_hyperslab(file_space, H5S_SELECT_SET, start, NULL, count, NULL);
#ifdef _PARTICLE
            H5Dwrite(partDataset_id[ispec], H5T_NATIVE_DOUBLE, partMemSpace, file_space, H5P_DEFAULT, &((*cuParticles)[ p ]->position(0)));
#endif
            H5Sclose(file_space);


        } // End for p

    } // End for ispec

}
