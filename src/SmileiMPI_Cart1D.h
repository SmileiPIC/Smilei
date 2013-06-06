
#ifndef SMILEIMPI_CART1D_H
#define SMILEIMPI_CART1D_H

#include "SmileiMPI.h"

#include <vector>
#include <string>

class Species;

class SmileiMPI_Cart1D : public SmileiMPI {
public:
	SmileiMPI_Cart1D( int* argc, char*** argv );
	virtual ~SmileiMPI_Cart1D();

	virtual void createTopology();
	virtual void exchangeParticles(Species* species, PicParams* params);

	inline int getProcCoord(int i) {return coords_[i];}

	inline bool isWester() {return (coords_[0]==0);}
	inline bool isEaster() {return (coords_[0]==dims_[0]-1);}

protected:
	MPI_Comm SMILEI_COMM_1D;

	int ndims_;
    int* dims_;

	// Cartesian ...
    int* coords_;
	int* periods_;
	int reorder_;

	int nbNeighbors_;
	int* neighbor_;

	virtual void sumFieldDual( Field* field );
	virtual void sumFieldPrim( Field* field );
	virtual void exchangeFieldDual( Field* field );
	virtual void exchangeFieldPrim( Field* field );
	virtual void writeFieldDual( Field* field, std::string name );
	virtual void writeFieldPrim( Field* field, std::string name );

};


/*
	if (n_particles!=0) {
		cout << "Size of particle : " << sizeof(*buff_send[0]) << endl;
		cout << "Diff between pointers" << endl;
		cout << &buff_send[0]->position(0)-(double*)(buff_send[0]) << endl;
		cout << &buff_send[0]->moments(0)-&buff_send[0]->position(0) << endl;
		//cout << &buff_send[0]->moments(0)-(double*)(buff_send[0]) << endl;
		cout << &buff_send[0]->moments(2)-&buff_send[0]->moments(0) << endl;
		cout << &buff_send[0]->chargeDensity()-&buff_send[0]->moments(2) << endl;
		//cout << &buff_send[0]->position(0)-&buff_send[0]->chargeDensity() << endl;

		cout << &buff_send[1]->position(0)-(double*)(buff_send[0]) << endl;
		cout << &buff_send[1]->position(0)-&buff_send[0]->position(0) << endl;
		cout << &buff_send[1]->moments(0)-&buff_send[1]->position(0) << endl;
		cout << &buff_send[1]->moments(2)-&buff_send[1]->moments(0) << endl;
		cout << &buff_send[1]->chargeDensity()-&buff_send[1]->moments(2) << endl;

		cout << "Diff between values" << endl;

		cout << buff_send[0]->moments(0) << " " << (double)(*((double*)(buff_send[0])+7)) << endl;
		cout << buff_send[0]->moments(1) << " " << (double)(*((double*)(buff_send[0])+8)) << endl;
		cout << buff_send[0]->moments(2) << " " << (double)(*((double*)(buff_send[0])+9)) << endl;
		cout << buff_send[0]->chargeDensity() << " " << (double)(*((double*)(buff_send[0])+10)) << endl;

		cout << buff_send[0]->position(0) << " " << (double)(*((double*)(buff_send[0])+6)) << endl;
		//cout << buff_send[0]->position(0) << " " << (double)(*((double*)(buff_send[0])+7)) << endl;
		//cout << buff_send[0]->position(0) << " " << (double)(*((double*)(buff_send[0])+8)) << endl;
	}
*/




#endif

