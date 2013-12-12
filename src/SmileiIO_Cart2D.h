/*
 * SmileIO_Cart2D.h
 *
 *  Created on: 3 juil. 2013
 *      Author: jderouil
 */

#ifndef SMILEIO_CART2D_H
#define SMILEIO_CART2D_H

#include "SmileiIO.h"


#include <string>
#include <vector>

class SmileiIO_Cart2D : public SmileiIO {
public:
	SmileiIO_Cart2D( PicParams* params, SmileiMPI* smpi );
	~SmileiIO_Cart2D();

	void open (  );
	void close(  );
	void createPattern( PicParams* params, SmileiMPI* smpi );

	void writeFieldsSingleFile( Field* field, hid_t file_id, int itime );
	void writeFieldsSingleFileTime( Field* field, hid_t group_id );

	void write( Field* field, std::string name );
	void write( Field* field, std::string name, double time );
	void writePerProcess( Field* field, std::string name, double time, int rank );


private:
	std::vector<unsigned int> istart;
	std::vector<unsigned int> bufsize;

	// [primDual][primDual]
	hid_t memspace_ [2][2];
	hid_t filespace_[2][2];

};





#endif /* SMILEIO_CART2D_H_ */
