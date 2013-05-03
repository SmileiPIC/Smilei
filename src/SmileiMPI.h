
#ifndef SMILEIMPI_H
#define SMILEIMPI_H

#include <mpi.h>

#include "PicParams.h"
#include <string>
#include <vector>

class PicParams;
class Species;
class Particle;

class ElectroMagn;
class Field;

class SmileiMPI {
public:
	SmileiMPI( int* argc, char*** argv );
	virtual ~SmileiMPI();

	inline bool isMaster() { return (smilei_rk==0); }
	inline void barrier() { MPI_Barrier( SMILEI_COMM_WORLD );}

	void bcast( PicParams& params );
	void bcast( std::string& val );
	void bcast( unsigned int &val );
	void bcast( double& val );
	void bcast( bool& val );
	void bcast( std::vector<unsigned int>& val );
	void bcast( std::vector<int>& val );
	void bcast( std::vector<double>& val );
	void bcast( SpeciesStructure& speciesStructure );
	void bcast( std::vector<SpeciesStructure>& vecSpeciesStructure );
	void bcast( LaserStructure& laserStructure );
	void bcast( std::vector<LaserStructure>& vecLaserStructure );
	void bcast_type_der( PicParams& params );

	virtual void createTopology() = 0;
	virtual void exchangeParticles(Species* species, PicParams* params) = 0;

	void sumRho( ElectroMagn* champs );
	void sumDensities( ElectroMagn* champs );

	virtual void sumField( Field* field ) = 0;

	void exchangeE( ElectroMagn* champs );
	void exchangeB( ElectroMagn* champs );

	virtual void exchangeField( Field* field ) = 0;

	inline int getRank() {return smilei_rk;}

 protected:
	MPI_Comm SMILEI_COMM_WORLD;

	int smilei_sz;
	int smilei_rk;

	std::vector<Particle*>* buff_send;
	std::vector<Particle*>* buff_recv;

	std::vector<unsigned int> oversize;
	std::vector<double> cell_length;

};

#endif

