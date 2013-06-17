
#ifndef SMILEIMPI_CART2D_H
#define SMILEIMPI_CART2D_H

#include "SmileiMPI.h"

#include <vector>
#include <string>

class Species;

class SmileiMPI_Cart2D : public SmileiMPI {
public:
	SmileiMPI_Cart2D( int* argc, char*** argv );
	SmileiMPI_Cart2D(SmileiMPI *smpi);
	virtual ~SmileiMPI_Cart2D();

	virtual void whoami() {std::cout << "SmileiMPI_Cart2D" << std::endl;}

	virtual void createTopology();
	virtual void exchangeParticles(Species* species, int ispec, PicParams* params);

	inline int getProcCoord(int i) {return coords_[i];}

	inline bool isWester() {return (coords_[0]==0);}
	inline bool isEaster() {return (coords_[0]==dims_[0]-1);}

	void createType( PicParams& params );

	virtual void exchangeField ( Field* field );
	virtual void sumField      ( Field* field );

protected:
	MPI_Comm SMILEI_COMM_2D;

	int ndims_;
    int* dims_;

	// Cartesian ...
    int* coords_;
	int* periods_;
	int reorder_;

	int nbNeighbors_;
	int neighbor_[3][2];	// 

	virtual void writeField    ( Field* field, std::string name );

	MPI_Datatype ntype_   [2][2][2];
	MPI_Datatype ntypeSum_[2][2][2];

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

