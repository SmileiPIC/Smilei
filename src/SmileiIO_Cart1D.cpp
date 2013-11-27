/*
 * SmileiIO_Cart1D.cpp
 *
 *  Created on: 3 juil. 2013
 *      Author: jderouil
 */

#include "SmileiIO_Cart1D.h"

#include "PicParams.h"
#include "SmileiMPI_Cart1D.h"
#include "Field1D.h"

using namespace std;

SmileiIO_Cart1D::SmileiIO_Cart1D( PicParams* params, SmileiMPI* smpi )
 : SmileiIO( params, smpi )
{
	open();
	createPattern(params,smpi);
}

SmileiIO_Cart1D::~SmileiIO_Cart1D()
{
	close();
}

void SmileiIO_Cart1D::open(  )
{
	MPI_Info info  = MPI_INFO_NULL;

    /*
     * Set up file access property list with parallel I/O access
     */
	hid_t plist_id = H5Pcreate(H5P_FILE_ACCESS);
    H5Pset_fapl_mpio(plist_id, MPI_COMM_WORLD, info);
   /*
    * Create a new file collectively and release property list identifier.
    */
    string fieldStr, dirStr, name;
	for (int itype=0 ; itype<4 ; itype++) {
		if      (itype==0) fieldStr = "e";
		else if (itype==1) fieldStr = "b";
		else if (itype==2) fieldStr = "j";
		else if (itype==3) {
			fieldStr = "r";
			name = "rho.h5";
			file_id_ [itype][0] = H5Fcreate( name.c_str(), H5F_ACC_TRUNC, H5P_DEFAULT, plist_id);
		}
		if (itype!=3) {
			for (int iDim=0 ; iDim<3 ; iDim++) {
				if      (iDim==0) dirStr = "x";
				else if (iDim==1) dirStr = "y";
				else if (iDim==2) dirStr = "z";
				name = fieldStr + dirStr + ".h5";
				file_id_ [itype][iDim] = H5Fcreate( name.c_str(), H5F_ACC_TRUNC, H5P_DEFAULT, plist_id);
			}
		}
	}
	H5Pclose(plist_id);

   //
   // Create property list for collective dataset write.
   //
	write_plist = H5Pcreate(H5P_DATASET_XFER);
	H5Pset_dxpl_mpio(write_plist, H5FD_MPIO_INDEPENDENT);

}

void SmileiIO_Cart1D::close(  )
{
//	for (int ix_isPrim=0 ; ix_isPrim<1 ; ix_isPrim++) {
//		for (int iy_isPrim=0 ; iy_isPrim<1 ; iy_isPrim++) {
//			H5Sclose( filespace_[1]);
//			H5Sclose( memspace_ [1]);
//		}
//	}

	for (int itype=0 ; itype<4 ; itype++) {
		if (itype!=3) {
			for (int iDim=0 ; iDim<3 ; iDim++) {
				H5Fclose( file_id_  [itype][iDim]);
			}
		}
		else
			H5Fclose( file_id_  [itype][0]);
	}


	H5Pclose( write_plist );
}

void SmileiIO_Cart1D::createPattern( PicParams* params, SmileiMPI* smpi )
{
	SmileiMPI_Cart1D* smpi1D =  static_cast<SmileiMPI_Cart1D*>(smpi);

	istart = smpi1D->oversize;
	bufsize.resize(params->nDim_field, 0);

	for (unsigned int i=0 ; i<params->nDim_field ; i++) {
		if (smpi1D->getProcCoord(i)!=0) istart[i]+=1;
		bufsize[i] = params->n_space[i] + 1;
	}

	int nx0 = bufsize[0];
	int ny0 = bufsize[1] ;

	// memspace [primDual]
	// fimespace[primDual]
	int nx;
	for (int ix_isPrim=0 ; ix_isPrim<2 ; ix_isPrim++) {
		nx = nx0 + ix_isPrim;

		istart = smpi1D->oversize;
		bufsize.resize(params->nDim_field, 0);

		for (unsigned int i=0 ; i<params->nDim_field ; i++) {
			if (smpi1D->getProcCoord(i)!=0) istart[i]+=1;
			bufsize[i] = params->n_space[i] + 1;
		}
		bufsize[0] += ix_isPrim;

		/*
		 * Create the dataspace for the dataset.
		 */
		hsize_t     chunk_dims[1];
		chunk_dims[0] = nx + 2*smpi1D->oversize[0] ;
		hid_t memspace  = H5Screate_simple(params->nDim_field, chunk_dims, NULL);

		hsize_t     offset[1];
		hsize_t     stride[1];
		hsize_t     count[1];

		offset[0] = istart[0];
		stride[0] = 1;

		if (smpi1D->number_of_procs[0] != 1) {
			if ( ix_isPrim == 0 ) {
				if (smpi1D->getProcCoord(0)!=0)
					bufsize[0]--;
			}
			else {
				if ( (smpi1D->coords_[0]!=0) && (smpi1D->coords_[0]!=smpi1D->number_of_procs[0]-1) )
					bufsize[0] -= 2;
				else
					bufsize[0] -= 1;
			}
		}
		count[0] = bufsize[0];
		H5Sselect_hyperslab(memspace, H5S_SELECT_SET, offset, stride, count, NULL);
		memspace_ [ ix_isPrim ] = memspace;

		// Each process defines dataset in memory and writes it to the hyperslab
		// in the file.
		//
		hsize_t     dimsf[1];
		dimsf[0] = params->n_space_global[0]+1+ix_isPrim;

		hid_t filespace = H5Screate_simple(params->nDim_field, dimsf, NULL);
		//
		// Select hyperslab in the file.
		//
		offset[0] = smpi->getCellStartingGlobalIndex(0)+istart[0];
		stride[0] = 1;
		count[0] = 1;
		hsize_t     block[1];
		block[0] = bufsize[0];
		H5Sselect_hyperslab(filespace, H5S_SELECT_SET, offset, stride, count, block);
		filespace_[ ix_isPrim ] = filespace;

	}

}

void SmileiIO_Cart1D::write( Field* field, string name )
{
	std::vector<unsigned int> isPrimal = field->isPrimal_;
	Field1D* f1D =  static_cast<Field1D*>(field);

	hid_t memspace  = memspace_ [ isPrimal[0] ];
	hid_t filespace = filespace_[ isPrimal[0] ];

	string typeStr = name.substr(0,1);
	string dirStr  = name.substr(1,1);
	int itype;
	if      ( typeStr.compare( (string)"e" ) == 0 ) itype = 0;
	else if ( typeStr.compare( (string)"b" ) == 0 ) itype = 1;
	else if ( typeStr.compare( (string)"j" ) == 0 ) itype = 2;
	else if ( typeStr.compare( (string)"r" ) == 0 ) itype = 3;

	int idir;
	if (itype != 3) {
		if      ( dirStr.compare( (string)"x" ) == 0 ) idir = 0;
		else if ( dirStr.compare( (string)"y" ) == 0 ) idir = 1;
		else if ( dirStr.compare( (string)"z" ) == 0 ) idir = 2;
	}
	else {
		idir = 0;
	}

	hid_t file_id   = file_id_  [ itype ][ idir ];

	hid_t plist_id = H5Pcreate(H5P_DATASET_CREATE);
//	chunk_dims[0] = bufsize[0];
//	H5Pset_chunk(plist_id, 2, chunk_dims); // Problem different dims for each process, may be the same for all ...
	hid_t dset_id = H5Dcreate(file_id, "Field", H5T_NATIVE_DOUBLE, filespace, H5P_DEFAULT, plist_id, H5P_DEFAULT);
	H5Pclose(plist_id);

	H5Dwrite( dset_id, H5T_NATIVE_DOUBLE, memspace, filespace, write_plist, &(f1D->data_[0]) );
	H5Dclose(dset_id);

}
