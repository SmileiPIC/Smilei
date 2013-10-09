/*
 * SmileiIO_Cart2D.cpp
 *
 *  Created on: 3 juil. 2013
 *      Author: jderouil
 */

#include "SmileiIO_Cart2D.h"

#include "PicParams.h"
#include "SmileiMPI_Cart2D.h"
#include "Field2D.h"

#include <sstream>

using namespace std;

SmileiIO_Cart2D::SmileiIO_Cart2D( PicParams* params, SmileiMPI* smpi )
 : SmileiIO( params, smpi )
{
	open();
	createPattern(params,smpi);
}

SmileiIO_Cart2D::~SmileiIO_Cart2D()
{
	close();
}

void SmileiIO_Cart2D::open(  )
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

void SmileiIO_Cart2D::close(  )
{
//	for (int ix_isPrim=0 ; ix_isPrim<1 ; ix_isPrim++) {
//		for (int iy_isPrim=0 ; iy_isPrim<1 ; iy_isPrim++) {
//			H5Sclose( filespace_[0][1]);
//			H5Sclose( memspace_ [0][1]);
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

void SmileiIO_Cart2D::createPattern( PicParams* params, SmileiMPI* smpi )
{
	SmileiMPI_Cart2D* smpi2D =  static_cast<SmileiMPI_Cart2D*>(smpi);

	istart = smpi2D->oversize;
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


//			if ( ( iy_isPrim == 0 ) && ( ix_isPrim == 1 ) ) {
//				cout << smpi->getRank() << " " << offset[0] << " " << offset[1] << " " << stride[0] << " " << stride[1] << " " << count[0] << " " << count[1] << endl;
//			}
			//
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

}

void SmileiIO_Cart2D::write( Field* field, string name, double time )
{
	MPI_Info info  = MPI_INFO_NULL;
	hid_t plist_id = H5Pcreate(H5P_FILE_ACCESS);
	H5Pset_fapl_mpio(plist_id, MPI_COMM_WORLD, info);
	ostringstream name_t;
	name_t.str(""); name_t << name << "_" << time << ".h5";
	hid_t file_t = H5Fcreate( name_t.str().c_str(), H5F_ACC_TRUNC, H5P_DEFAULT, plist_id);
	H5Pclose(plist_id);

	std::vector<unsigned int> isPrimal = field->isPrimal_;
	Field2D* f2D =  static_cast<Field2D*>(field);

	hid_t memspace  = memspace_ [ isPrimal[0] ][ isPrimal[1] ];
	hid_t filespace = filespace_[ isPrimal[0] ][ isPrimal[1] ];

	plist_id = H5Pcreate(H5P_DATASET_CREATE);
	hid_t dset_id = H5Dcreate(file_t, "Field", H5T_NATIVE_DOUBLE, filespace, H5P_DEFAULT, plist_id, H5P_DEFAULT);
	H5Pclose(plist_id);
	H5Dwrite( dset_id, H5T_NATIVE_DOUBLE, memspace, filespace, write_plist, &(f2D->data_[0][0]) );
	H5Dclose( dset_id );


	H5Fclose( file_t );

}

void SmileiIO_Cart2D::write( Field* field, string name )
{
	std::vector<unsigned int> isPrimal = field->isPrimal_;
	Field2D* f2D =  static_cast<Field2D*>(field);

	hid_t memspace  = memspace_ [ isPrimal[0] ][ isPrimal[1] ];
	hid_t filespace = filespace_[ isPrimal[0] ][ isPrimal[1] ];

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
//	chunk_dims[1] = bufsize[1];
//	H5Pset_chunk(plist_id, 2, chunk_dims); // Problem different dims for each process, may be the same for all ...
	hid_t dset_id = H5Dcreate(file_id, "Field", H5T_NATIVE_DOUBLE, filespace, H5P_DEFAULT, plist_id, H5P_DEFAULT);
	H5Pclose(plist_id);

	H5Dwrite( dset_id, H5T_NATIVE_DOUBLE, memspace, filespace, write_plist, &(f2D->data_[0][0]) );
	H5Dclose(dset_id);


}


void SmileiIO_Cart2D::writePerProcess( Field* field, string name, double time, int rank )
{
	ostringstream iname;

	iname.str(""); iname << name << "_" << rank << ".h5";

	Field2D* f2D =  static_cast<Field2D*>(field);
	hid_t       file_id;
	herr_t      status;

	file_id = H5Fcreate(iname.str().c_str(), H5F_ACC_TRUNC, H5P_DEFAULT, H5P_DEFAULT);

	hid_t       dataset_id, dataspace_id;
	hsize_t     dims[2];
	dims[0] = field->dims_[0];
	dims[1] = field->dims_[1];
	dataspace_id = H5Screate_simple(2, dims, NULL);
	dataset_id = H5Dcreate2(file_id, "Field", H5T_NATIVE_DOUBLE, dataspace_id, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);

	status = H5Dwrite(dataset_id, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, &(f2D->data_[0][0]));

	status = H5Dclose(dataset_id);
	status = H5Sclose(dataspace_id);

	status = H5Fclose(file_id);

}
