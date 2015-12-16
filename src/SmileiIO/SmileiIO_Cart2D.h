/*
 * SmileIO_Cart2D.h
 *
 *  Created on: 3 juil. 2013
 */
#ifndef SMILEIO_CART2D_H
#define SMILEIO_CART2D_H

#include <string>
#include <vector>

#include "SmileiIO.h"

//  --------------------------------------------------------------------------------------------------------------------
//! Class SmileiIO_Cart2D
//  --------------------------------------------------------------------------------------------------------------------
class SmileiIO_Cart2D : public SmileiIO {
public:
    //! Create // HDF5 environment
    SmileiIO_Cart2D( Params& params, Diagnostic& diag, SmileiMPI* smpi );
    //! Destructor for SmileiIO
    ~SmileiIO_Cart2D();

    //! Build memory and file space for // HDF5 write/read
    void createPattern( Params& params, SmileiMPI* smpi );

    //! Basic write current field in specified group of the global file
    void writeFieldsSingleFileTime( Field* field, hid_t group_id );
    void writeOneFieldSingleFileTime( Field* field, hid_t group_id );

    //! Basic write field on its own file (debug)
    void write( Field* field );

private:
    //! memory space for // HDF5 write/read
    //! Size = 2 x 2 : 0 if prim, 1 if dual per direction
    hid_t memspace_ [2][2];
    //! file space for // HDF5 write/read
    //! Size = 2 x 2 : 0 if prim, 1 if dual per direction
    hid_t filespace_[2][2];

    //! \todo Define chunk size of output for interpolated output
    //hsize_t chunk_dims[2];

};

#endif /* SMILEIO_CART2D_H_ */
