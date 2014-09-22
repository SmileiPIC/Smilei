#include "SmileiMPI.h"

#include <cmath>
#include <cstring>

#include <iostream>
#include <sstream>

#include "PicParams.h"
#include "DiagParams.h"
#include "Tools.h"

#include "ElectroMagn.h"
#include "Field.h"

#include "Species.h"

using namespace std;

SmileiMPI::SmileiMPI( int* argc, char*** argv )
{    
    int mpi_provided;

    MPI_Init_thread( argc, argv, MPI_THREAD_FUNNELED, &mpi_provided );
    if (mpi_provided == MPI_THREAD_SINGLE){
        MESSAGE("openMP not supported");
    }

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

/*
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

    for(int k=0; k<vecSize; k++) {
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
    bcast( speciesStructure.bc_part_type_long );
    bcast( speciesStructure.bc_part_type_trans );
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
    bcast( laserStructure.y_profile );
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
*/
void SmileiMPI::sumRho( ElectroMagn* EMfields )
{
    sumField( EMfields->rho_ );

}

void SmileiMPI::sumRhoJ( ElectroMagn* EMfields )
{
    // sum total charge density and currents
    sumField( EMfields->rho_ );
    sumField( EMfields->Jx_ );
    sumField( EMfields->Jy_ );
    sumField( EMfields->Jz_ );

    // sum density and currents for all species
    for (unsigned int ispec=0; ispec<EMfields->n_species; ispec++) {
        sumField( EMfields->rho_s[ispec] );
        sumField( EMfields->Jx_s[ispec] );
        sumField( EMfields->Jy_s[ispec] );
        sumField( EMfields->Jz_s[ispec] );
    }

}

void SmileiMPI::exchangeE( ElectroMagn* EMfields )
{
    exchangeField( EMfields->Ex_ );
    exchangeField( EMfields->Ey_ );
    exchangeField( EMfields->Ez_ );

}

void SmileiMPI::exchangeB( ElectroMagn* EMfields )
{
    exchangeField( EMfields->Bx_ );
    exchangeField( EMfields->By_ );
    exchangeField( EMfields->Bz_ );

}
