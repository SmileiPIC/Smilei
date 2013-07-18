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
	for (int itype=0 ; itype<2 ; itype++) {
		if      (itype==0) fieldStr = "e";
		else if (itype==1) fieldStr = "b";
		for (int iDim=0 ; iDim<3 ; iDim++) {
			if      (iDim==0) dirStr = "x";
			else if (iDim==1) dirStr = "y";
			else if (iDim==2) dirStr = "z";
			name = fieldStr + dirStr+".h5";
			file_id_ [itype][iDim] = H5Fcreate( name.c_str(), H5F_ACC_TRUNC, H5P_DEFAULT, plist_id);
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

	for (int itype=0 ; itype<2 ; itype++) {
		for (int iDim=0 ; iDim<3 ; iDim++) {
			H5Fclose( file_id_  [itype][iDim]);
		}
	}


	H5Pclose( write_plist );
}

void SmileiIO_Cart1D::createPattern( PicParams* params, SmileiMPI* smpi )
{

}

void SmileiIO_Cart1D::write( Field* field, string name, double time )
{

}
