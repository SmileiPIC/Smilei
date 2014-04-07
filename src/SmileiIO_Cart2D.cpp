/*
 * SmileiIO_Cart2D.cpp
 *
 *  Created on: 3 juil. 2013
 *      Author: jderouil
 */

#include "SmileiIO_Cart2D.h"

#include <sstream>

#include "PicParams.h"
#include "SmileiMPI_Cart2D.h"
#include "Field2D.h"

using namespace std;

SmileiIO_Cart2D::SmileiIO_Cart2D( PicParams* params, SmileiMPI* smpi )
    : SmileiIO( params, smpi )
{
    createPattern(params,smpi);
}

SmileiIO_Cart2D::~SmileiIO_Cart2D()
{
}

void SmileiIO_Cart2D::createPattern( PicParams* params, SmileiMPI* smpi )
{
    SmileiMPI_Cart2D* smpi2D =  static_cast<SmileiMPI_Cart2D*>(smpi);

    std::vector<unsigned int> istart;
    istart = smpi2D->oversize;
    std::vector<unsigned int> bufsize;
    bufsize.resize(params->nDim_field, 0);

    for (unsigned int i=0 ; i<params->nDim_field ; i++) {
        if (smpi2D->getProcCoord(i)!=0) istart[i]+=1;
        bufsize[i] = params->n_space[i] + 1;
    }

    int nx0 = bufsize[0];
    int ny0 = bufsize[1] ;

    // memspace [primDual][primDual]
    // fimespace[primDual][primDual]
    int nx, ny;
    for (int ix_isPrim=0 ; ix_isPrim<2 ; ix_isPrim++) {
        nx = nx0 + ix_isPrim;
        for (int iy_isPrim=0 ; iy_isPrim<2 ; iy_isPrim++) {
            ny = ny0 + iy_isPrim;

            istart = smpi2D->oversize;
            bufsize.resize(params->nDim_field, 0);

            for (unsigned int i=0 ; i<params->nDim_field ; i++) {
                if (smpi2D->getProcCoord(i)!=0) istart[i]+=1;
                bufsize[i] = params->n_space[i] + 1;
            }
            bufsize[0] += ix_isPrim;
            bufsize[1] += iy_isPrim;


            /*
             * Create the dataspace for the dataset.
             */
            hsize_t     chunk_dims[2];
            chunk_dims[0] = nx + 2*smpi2D->oversize[0] ;
            chunk_dims[1] = ny + 2*smpi2D->oversize[1] ;
            hid_t memspace  = H5Screate_simple(params->nDim_field, chunk_dims, NULL);

            hsize_t     offset[2];
            hsize_t     stride[2];
            hsize_t     count[2];

            offset[0] = istart[0];
            offset[1] = istart[1];
            stride[0] = 1;
            stride[1] = 1;

            if (smpi2D->number_of_procs[0] != 1) {
                if ( ix_isPrim == 0 ) {
                    if (smpi2D->getProcCoord(0)!=0)
                        bufsize[0]--;
                }
                else {
                    if ( (smpi2D->coords_[0]!=0) && (smpi2D->coords_[0]!=smpi2D->number_of_procs[0]-1) )
                        bufsize[0] -= 2;
                    else
                        bufsize[0] -= 1;
                }
            }
            if (smpi2D->number_of_procs[1] != 1) {
                if ( iy_isPrim == 0 ) {
                    if (smpi2D->getProcCoord(1)!=0)
                        bufsize[1]--;
                }
                else {
                    if ( (smpi2D->coords_[1]!=0) && (smpi2D->coords_[1]!=smpi2D->number_of_procs[1]-1) )
                        bufsize[1] -= 2;
                    else
                        bufsize[1] -= 1;
                }
            }
            count[0] = bufsize[0];
            count[1] = bufsize[1];
            H5Sselect_hyperslab(memspace, H5S_SELECT_SET, offset, stride, count, NULL);
            memspace_ [ ix_isPrim ][ iy_isPrim ] = memspace;


	    //if ( ( iy_isPrim == 0 ) && ( ix_isPrim == 1 ) ) {
	    //    cout << smpi->getRank() << " " << offset[0] << " " << offset[1] << " " << stride[0] << " " << stride[1] << " " << count[0] << " " << count[1] << endl;
	    //}

            // Each process defines dataset in memory and writes it to the hyperslab
            // in the file.
            //
            hsize_t     dimsf[2];
            dimsf[0] = params->n_space_global[0]+1+ix_isPrim;
            dimsf[1] = params->n_space_global[1]+1+iy_isPrim;

            hid_t filespace = H5Screate_simple(params->nDim_field, dimsf, NULL);
            //
            // Select hyperslab in the file.
            //
            offset[0] = smpi->getCellStartingGlobalIndex(0)+istart[0];
            offset[1] = smpi->getCellStartingGlobalIndex(1)+istart[1];
            stride[0] = 1;
            stride[1] = 1;
            count[0] = 1;
            count[1] = 1;
            hsize_t     block[2];
            block[0] = bufsize[0];
            block[1] = bufsize[1];
            H5Sselect_hyperslab(filespace, H5S_SELECT_SET, offset, stride, count, block);
            filespace_[ ix_isPrim ][ iy_isPrim ] = filespace;

        }
    }

} // END createPattern


// ---------------------------------------------------------------------------------------------------------------------
// For time iteration "itime", write current field in the time step dataset of the global file
// ---------------------------------------------------------------------------------------------------------------------
void SmileiIO_Cart2D::writeFieldsSingleFileTime( Field* field, hid_t group_id )
{
    std::vector<unsigned int> isDual = field->isDual_;
    Field2D* f2D =  static_cast<Field2D*>(field);

    hid_t memspace  = memspace_ [ isDual[0] ][ isDual[1] ];
    hid_t filespace = filespace_[ isDual[0] ][ isDual[1] ];

    hid_t plist_id = H5Pcreate(H5P_DATASET_CREATE);
    //H5Pset_chunk(plist_id, 2, chunk_dims); // Problem different dims for each process
    //hid_t dset_id = H5Dcreate(file_id, (field->name).c_str(), H5T_NATIVE_DOUBLE, filespace, H5P_DEFAULT, plist_id, H5P_DEFAULT);
    hid_t dset_id = H5Dcreate2(group_id, (field->name).c_str(), H5T_NATIVE_DOUBLE, filespace, H5P_DEFAULT, plist_id, H5P_DEFAULT);

    H5Pclose(plist_id);

    H5Dwrite( dset_id, H5T_NATIVE_DOUBLE, memspace, filespace, write_plist, &(f2D->data_2D[0][0]) );
    H5Dclose(dset_id);


} // END writeFieldsSingleFileTime


//! this method writes a field on an hdf5 file should be used just for debug
void SmileiIO_Cart2D::write( Field* field )
{
    std::vector<unsigned int> isDual = field->isDual_;
    Field2D* f2D =  static_cast<Field2D*>(field);
    string name = f2D->name+".h5";

    hid_t memspace  = memspace_ [ isDual[0] ][ isDual[1] ];
    hid_t filespace = filespace_[ isDual[0] ][ isDual[1] ];

    MPI_Info info  = MPI_INFO_NULL;
    hid_t plist_id = H5Pcreate(H5P_FILE_ACCESS);
    H5Pset_fapl_mpio(plist_id, MPI_COMM_WORLD, info);
    hid_t file_id = H5Fcreate( name.c_str(), H5F_ACC_TRUNC, H5P_DEFAULT, plist_id);
    H5Pclose(plist_id);

    plist_id = H5Pcreate(H5P_DATASET_CREATE);
    // H5Pset_chunk(plist_id, 2, chunk_dims); // Problem different dims for each process, may be the same for all ...
    hid_t dset_id = H5Dcreate(file_id, "Field", H5T_NATIVE_DOUBLE, filespace, H5P_DEFAULT, plist_id, H5P_DEFAULT);
    H5Pclose(plist_id);

    H5Dwrite( dset_id, H5T_NATIVE_DOUBLE, memspace, filespace, write_plist, &(f2D->data_2D[0][0]) );
    H5Dclose(dset_id);

    H5Fclose( file_id );

} // END write

