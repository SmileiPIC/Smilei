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

    // To be used when python tools updated
    void writeAllFields( ElectroMagn* EMfields, int itime );
    virtual void writeFieldsSingleFile( Field* field, hid_t file_id, int itime ) = 0;

    // To be used when python tools updated
    void writeAllFieldsSingleFileTime( ElectroMagn* EMfields, int itime );
    virtual void writeFieldsSingleFileTime( Field* field, hid_t group_id ) = 0;
    hid_t global_file_id_;

    void writePlasma( std::vector<Species*> vecSpecies, double time, SmileiMPI* smpi );


    // [2][3] : [2] -> 0 = E, 1 = B
    // [2][3] : [3] -> x, y, z
    hid_t file_id_  [4][3];

    hid_t write_plist;

    hid_t  partFile_id;
    unsigned int nDatasetSpecies;
    hid_t* partDataset_id;  /* identifiers */
    hid_t partMemSpace;
    int particleSize;

    //! for debugging
    virtual void write( Field* field ) = 0;

private:
};

#endif /* SMILEI_OUTPUT_H_ */
