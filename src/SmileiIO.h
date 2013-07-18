/*
 * SmileiIO.h
 *
 *  Created on: 3 juil. 2013
 *      Author: jderouil
 */

#ifndef SMILEIIO_H
#define SMILEIIO_H

#include <string>
#include <vector>
#include <hdf5.h>

class PicParams;
class SmileiMPI;
class ElectroMagn;
class Field;
class Species;

class SmileiIO {
public:
	SmileiIO( PicParams* params, SmileiMPI* smpi );
	virtual ~SmileiIO();

	virtual void write( Field* field, std::string name, double time ) = 0;

	void writeFields( ElectroMagn* EMfields, double time );
	void writePlasma( std::vector<Species*> vecSpecies, double time, SmileiMPI* smpi );


	// [2][3] : [2] -> 0 = E, 1 = B
	// [2][3] : [3] -> x, y, z
	hid_t file_id_  [2][3];

	hid_t write_plist;

	hid_t  partFile_id;
	int nDatasetSpecies;
	hid_t* partDataset_id;  /* identifiers */
	hid_t partMemSpace;


private:
};

#endif /* SMILEI_OUTPUT_H_ */
