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

	void open (  );
	void close(  );
	void createPattern( PicParams* params, SmileiMPI* smpi );

	void writeFieldsSingleFile( Field* field, hid_t file_id, int itime );
	void writeFieldsSingleFileTime( Field* field, hid_t group_id );

	void write( Field* field, std::string name );
	void write( Field* field, std::string name, double time ) {MESSAGE("not implemented");}
	void writePerProcess( Field* field, std::string name, double time, int rank ){MESSAGE("not implemented");}


private:
	std::vector<unsigned int> istart;
	std::vector<unsigned int> bufsize;

	// [primDual][primDual]
	hid_t memspace_ [2];
	hid_t filespace_[2];

};





#endif /* SMILEIO_CART1D_H_ */
