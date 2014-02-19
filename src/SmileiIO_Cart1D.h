/*
 * SmileIO_Cart1D.h
 *
 *  Created on: 3 juil. 2013
 *      Author: jderouil
 */

#ifndef SMILEIO_CART1D_H
#define SMILEIO_CART1D_H

#include "SmileiIO.h"
#include "Tools.h"

#include <string>
#include <vector>

class SmileiIO_Cart1D : public SmileiIO {
public:
    SmileiIO_Cart1D( PicParams* params, SmileiMPI* smpi );
    ~SmileiIO_Cart1D();

    void createPattern( PicParams* params, SmileiMPI* smpi );

    void writeFieldsSingleFile( Field* field, hid_t file_id, int itime );
    void writeFieldsSingleFileTime( Field* field, hid_t group_id );

    //! for debugging
    void write( Field* field );

private:
    std::vector<unsigned int> istart;
    std::vector<unsigned int> bufsize;

    // [primDual][primDual]
    hid_t memspace_ [2];
    hid_t filespace_[2];

};





#endif /* SMILEIO_CART1D_H_ */
