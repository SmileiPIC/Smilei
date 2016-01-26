/*
 * SmileIO_Cart1D.h
 *
 *  Created on: 3 juil. 2013
 */

#ifndef SMILEIO_CART1D_H
#define SMILEIO_CART1D_H

#include <string>
#include <vector>

#include "SmileiIO.h"
#include "Tools.h"

//  --------------------------------------------------------------------------------------------------------------------
//! Class SmileiIO_Cart1D
//  --------------------------------------------------------------------------------------------------------------------
class SmileiIO_Cart1D : public SmileiIO {
public:
    //! Create // HDF5 environment
    SmileiIO_Cart1D( Params& params, Diagnostic *diag, Patch* patch );
    //! Destructor for SmileiIO
    ~SmileiIO_Cart1D();

    //! Build memory and file space for // HDF5 write/read
    void createPattern( Params& params, Patch* patch );
    void updatePattern( Params& params, Patch* patch );

    //! Basic write current field in specified group of the global file
    void writeFieldsSingleFileTime( Field* field, hid_t group_id );
    void writeOneFieldSingleFileTime( Field* field, hid_t group_id ) {;}

    //! Basic write field on its own file (debug)
    void write( Field* field );

private:
    //! memory space for // HDF5 write/read
    //! Size = 2 : 0 if prim, 1 if dual
    hid_t memspace_ [2];
    //! file space for // HDF5 write/read
    //! Size = 2 : 0 if prim, 1 if dual
    hid_t filespace_[2];

    //! \todo Define chunk size of output for interpolated output
    //hsize_t chunk_dims[1];

};

#endif /* SMILEIO_CART1D_H */
