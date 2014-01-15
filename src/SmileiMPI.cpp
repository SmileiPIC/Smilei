#include "SmileiMPI.h"

#include "PicParams.h"
#include "DiagParams.h"
#include "Tools.h"

#include "ElectroMagn.h"
#include "Field.h"

#include "Species.h"

#include <iostream>
#include <sstream>
#include <cmath>
using namespace std;

SmileiMPI::SmileiMPI( int* argc, char*** argv )
{
	MPI_Init( argc, argv );
	SMILEI_COMM_WORLD = MPI_COMM_WORLD;
	MPI_Comm_size( SMILEI_COMM_WORLD, &smilei_sz );
	MPI_Comm_rank( SMILEI_COMM_WORLD, &smilei_rk );
    
}

SmileiMPI::SmileiMPI( SmileiMPI *smpi )
{
	SMILEI_COMM_WORLD = smpi->SMILEI_COMM_WORLD;
	MPI_Comm_size( SMILEI_COMM_WORLD, &smilei_sz );
	MPI_Comm_rank( SMILEI_COMM_WORLD, &smilei_rk );
    
	oversize = smpi->oversize;
	cell_starting_global_index = smpi->cell_starting_global_index;
	min_local = smpi->min_local;
	max_local = smpi->max_local;
    
	n_space_global = smpi->n_space_global;
    
}

SmileiMPI::~SmileiMPI()
{
	int status = 0;
	MPI_Finalized( &status );
	if (!status) MPI_Finalize();
    
}

void SmileiMPI::bcast( InputData& idata )
{
	DEBUG(10,"broadcast namelist");
	bcast(idata.namelist);
}

void SmileiMPI::bcast( string& val )
{
	int charSize=0;
	if (isMaster()) charSize = val.size()+1;
	MPI_Bcast(&charSize, 1, MPI_INT, 0, SMILEI_COMM_WORLD);
	
	char tmp[charSize];
	strcpy(tmp, val.c_str());
	MPI_Bcast(&tmp, charSize, MPI_CHAR, 0, SMILEI_COMM_WORLD);
    //	stringstream sstream;
    //	sstream << tmp;
    //	sstream >> val;
	
	if (!isMaster()) val=tmp;
	
}

void SmileiMPI::init( PicParams& params )
{
	oversize.resize(params.nDim_field, 0);
	cell_starting_global_index.resize(params.nDim_field, 0);
	min_local.resize(params.nDim_field, 0.);
	max_local.resize(params.nDim_field, 0.);
	n_space_global.resize(params.nDim_field, 0);
}


void SmileiMPI::bcast( short& val )
{
	MPI_Bcast( &val, 1, MPI_SHORT, 0, SMILEI_COMM_WORLD);
}

void SmileiMPI::bcast( unsigned int& val )
{
	MPI_Bcast( &val, 1, MPI_UNSIGNED, 0, SMILEI_COMM_WORLD);
}

void SmileiMPI::bcast( double& val )
{
	MPI_Bcast( &val, 1, MPI_DOUBLE, 0, SMILEI_COMM_WORLD);
}

void SmileiMPI::bcast( bool& val )
{
	MPI_Bcast( &val, 1, MPI_LOGICAL, 0, SMILEI_COMM_WORLD);
}

void SmileiMPI::bcast( vector<unsigned int>& val )
{
	int vecSize;
	if (isMaster()) vecSize = val.size();
	MPI_Bcast( &vecSize, 1, MPI_INT, 0, SMILEI_COMM_WORLD);
	if (!isMaster()) val.resize( vecSize );
    
	MPI_Bcast( &(val[0]), vecSize, MPI_UNSIGNED, 0, SMILEI_COMM_WORLD);
}

void SmileiMPI::bcast( vector<int>& val )
{
	int vecSize;
	if (isMaster()) vecSize = val.size();
	MPI_Bcast( &vecSize, 1, MPI_INT, 0, SMILEI_COMM_WORLD);
	if (!isMaster()) val.resize( vecSize );
    
	MPI_Bcast( &(val[0]), vecSize, MPI_INT, 0, SMILEI_COMM_WORLD);
}

void SmileiMPI::bcast( vector<double>& val )
{
	int vecSize=0;
	if (isMaster()) vecSize = val.size();
	MPI_Bcast( &vecSize, 1, MPI_INT, 0, SMILEI_COMM_WORLD);
	if (!isMaster()) val.resize( vecSize );
    
	MPI_Bcast( &(val[0]), vecSize, MPI_DOUBLE, 0, SMILEI_COMM_WORLD);
}

void SmileiMPI::bcast( vector<vector<double> >& val )
{
    int vecSize=0;
    if (isMaster()) vecSize = val.size();
    MPI_Bcast( &vecSize, 1, MPI_INT, 0, SMILEI_COMM_WORLD);
    if (!isMaster()) val.resize( vecSize );
    
    for(unsigned int k=0;k<vecSize; k++) {
        bcast(val[k]);
    }
    
}

void SmileiMPI::bcast( SpeciesStructure& speciesStructure )
{
	bcast( speciesStructure.species_type );
	bcast( speciesStructure.atomic_number );
	bcast( speciesStructure.initialization_type );
	bcast( speciesStructure.n_part_per_cell );
	bcast( speciesStructure.c_part_max );
	bcast( speciesStructure.mass );
	bcast( speciesStructure.charge );
	bcast( speciesStructure.density );
	bcast( speciesStructure.mean_velocity ); // must be params.nDim_field
	bcast( speciesStructure.temperature );
	bcast( speciesStructure.dynamics_type );
	bcast( speciesStructure.bc_part_type );
	bcast( speciesStructure.time_frozen );
	bcast( speciesStructure.radiating );
	bcast( speciesStructure.ionization_model );
	cout << "ioni ... " << speciesStructure.ionization_model << endl;
}

void SmileiMPI::bcast( vector<SpeciesStructure>& val )
{
	int vecSize;
	if (isMaster()) vecSize = val.size();
	MPI_Bcast( &vecSize, 1, MPI_INT, 0, SMILEI_COMM_WORLD);
	if (!isMaster()) val.resize( vecSize );
    
	for (int i=0 ; i<vecSize ; i++) bcast( val[i] );
}

void SmileiMPI::bcast( LaserStructure& laserStructure )
{
	bcast( laserStructure.a0 );
	bcast( laserStructure.angle );
	bcast( laserStructure.delta );
	bcast( laserStructure.time_profile );
	bcast( laserStructure.int_params );
	bcast( laserStructure.double_params );
}

void SmileiMPI::bcast( vector<LaserStructure>& val )
{
	int vecSize;
	if (isMaster()) vecSize = val.size();
	MPI_Bcast( &vecSize, 1, MPI_INT, 0, SMILEI_COMM_WORLD);
	if (!isMaster()) val.resize( vecSize );
    
	for (int i=0 ; i<vecSize ; i++) bcast( val[i] );
}

void SmileiMPI::bcast_type_der( PicParams& params )
{
	//! \todo{Using data structures}
}

void SmileiMPI::sumRho( ElectroMagn* champs )
{
    sumField( champs->rho_ );
    
}

void SmileiMPI::sumRhoJ( ElectroMagn* champs )
{
    // sum total charge density and currents
    sumField( champs->rho_ );
    sumField( champs->Jx_ );
    sumField( champs->Jy_ );
    sumField( champs->Jz_ );
    
    // sum density and currents for all species
    for (unsigned int ispec; ispec<champs->n_species; ispec++){
        sumField( champs->rho_s[ispec] );
        sumField( champs->Jx_s[ispec] );
        sumField( champs->Jy_s[ispec] );
        sumField( champs->Jz_s[ispec] );
    }
    
}

void SmileiMPI::exchangeE( ElectroMagn* champs )
{
    exchangeField( champs->Ex_ );
    exchangeField( champs->Ey_ );
    exchangeField( champs->Ez_ );
    
}

void SmileiMPI::exchangeB( ElectroMagn* champs )
{
    exchangeField( champs->Bx_ );
    exchangeField( champs->By_ );
    exchangeField( champs->Bz_ );
    
}

void SmileiMPI::solvePoissonPara( ElectroMagn* champs )
{
	for ( int i_rk = 0 ; i_rk < smilei_sz ; i_rk++ ) {
		if (i_rk==smilei_rk)
			champs->solvePoisson(this);
        
		barrier();
		exchangeField( champs->Ex_ );
	}
    
} // END solvePoissonPara


void SmileiMPI::writeFields( ElectroMagn* champs )
{
	writeField( champs->Ex_, "fex" );
	writeField( champs->Ey_, "fey" );
	writeField( champs->Ez_, "fez" );
	writeField( champs->Bx_, "fbx" );
	writeField( champs->By_, "fby" );
	writeField( champs->Bz_, "fbz" );
	writeField( champs->Jx_, "fjx" );
	writeField( champs->Jy_, "fjy" );
	writeField( champs->Jz_, "fjz" );
	writeField( champs->rho_, "rho" );
    
} // END writeFields

void SmileiMPI::writePlasma( vector<Species*> vecSpecies, string name )
{
	ofstream ofile;
	int n_species = vecSpecies.size();
    
	for (int ispec=0 ; ispec<n_species ; ispec++) {
        
		for ( int i_rk = 0 ; i_rk < smilei_sz ; i_rk++ ) {
			if (i_rk==smilei_rk) {
                if ((smilei_rk==0)&&(ispec==0)) ofile.open(name.c_str(), ios::out);
                else                            ofile.open(name.c_str(), ios::app);
                
                vecSpecies[ispec]->dump(ofile);
                ofile.close();
			}
			barrier();
		}
        
		ofile << endl;
	}
    
} // END writePlasma

